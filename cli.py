"""
kirsolve - resolucao de alelos KIR ambiguos a partir de saidas PING + kir-mapper.

Uso tipico:

    python -m kirsolve run dados.xlsx -o resultados/
    python -m kirsolve inspect dados.xlsx
    python -m kirsolve run dados.xlsx -o out/ --priors freq_brasil.tsv --prior-weight 5
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from .loaders import DEFAULT_SAMPLE_REGEX, load_sources, load_workbook
from .nomenclature import Allele
from .priors import (
    catalog_as_prior,
    fetch_ipd_catalog,
    load_prior_file,
    validate_against_catalog,
)
from .report import write_outputs
from .resolve import (
    consolidate,
    cross_tool,
    resolution_recommendation,
    resolve_cohort,
    ubiquity_report,
)


def cmd_inspect(args) -> int:
    """Mostra como cada arquivo ou diretorio foi interpretado."""
    from pathlib import Path as _P

    entradas = args.input if isinstance(args.input, list) else [args.input]
    for e in entradas:
        alvo = _P(e)
        if alvo.is_dir():
            achados = sorted(
                f for f in alvo.rglob("*")
                if f.suffix.lower() in {".csv", ".tsv", ".txt", ".xlsx", ".xls"}
            )
            print(f"Diretorio: {e}  ({len(achados)} arquivos legiveis)")
            for f in achados[:15]:
                print(f"   {f.relative_to(alvo)}")
            if len(achados) > 15:
                print(f"   ... e mais {len(achados)-15}")
        elif alvo.suffix.lower() in {".xlsx", ".xls"}:
            print(f"Planilha: {e}")
            print(f"   abas: {pd.ExcelFile(alvo).sheet_names}")
        else:
            print(f"Arquivo: {e}")
    print()

    df = load_sources(entradas, tool_map=_parse_sheet_tools(args.sheet_tool),
                      sample_regex=args.sample_regex)
    if df.empty:
        print("Nenhuma chamada reconhecida. Verifique --sheet-tool e --sample-regex.")
        return 1
    print(df.groupby(["source_sheet", "tool"]).size().to_string(), "\n")
    print("Amostras:", sorted(df["sample"].unique()), "\n")
    print("Loci:", sorted(df["locus"].unique()), "\n")
    warn = df[df["warnings"] != ""]
    if len(warn):
        print(f"{len(warn)} linhas com avisos de parsing (primeiras 10):")
        print(warn[["sample", "tool", "locus", "warnings"]].head(10).to_string(index=False))
    return 0


def cmd_run(args) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[1/6] lendo planilha...")
    raw = load_sources(args.input, tool_map=_parse_sheet_tools(args.sheet_tool),
                       sample_regex=args.sample_regex)
    if raw.empty:
        print("ERRO: nenhuma chamada reconhecida.", file=sys.stderr)
        return 1
    print(f"      {len(raw)} registros | {raw['sample'].nunique()} amostras "
          f"| {raw['locus'].nunique()} loci")

    print("[2/6] consolidando evidencias por ferramenta...")
    consolidated = consolidate(raw)

    print("[3/6] cruzando PING x kir-mapper...")
    cross = cross_tool(consolidated)
    print("      " + cross["status"].value_counts().to_dict().__str__())

    print("[4/6] montando priors...")
    priors = {}
    if args.priors:
        priors = load_prior_file(args.priors)
        print(f"      prior de arquivo: {sum(len(v) for v in priors.values())} alelos")
    elif args.use_ipd:
        observed = {}
        for _, r in cross.iterrows():
            observed.setdefault(r["locus"], set()).update(
                a.bare() for g in r["candidates"] for a in g if a.is_typed
            )
        catalog = fetch_ipd_catalog(outdir / "ipd_kir_allelelist.txt")
        priors = catalog_as_prior(catalog, observed)
        if catalog:
            observed_all = {k: v for k, v in observed.items()}
            unknown = validate_against_catalog(observed_all, catalog)
            if len(unknown):
                unknown.to_csv(outdir / "alelos_fora_do_catalogo.tsv", sep="\t", index=False)
                print(f"      {len(unknown)} alelos fora do catalogo IPD-KIR")

    print("[5/6] rodando EM e resolvendo...")
    calls, freqs, diags = resolve_cohort(
        cross, priors=priors, prior_weight=args.prior_weight,
        min_samples_for_em=args.min_samples_em,
    )
    ubiq = ubiquity_report(cross)
    resol = resolution_recommendation(cross)

    n_em = max((d["n_samples"] for d in diags.values()), default=0)
    if n_em < args.min_samples_em:
        print(f"      AVISO: no maximo {n_em} amostras informativas por locus. "
              f"O EM nao e confiavel abaixo de {args.min_samples_em}; as chamadas "
              f"'provavel_EM' foram suprimidas e so valem as resolucoes deterministicas.")

    print("[6/6] gravando relatorios...")
    files = write_outputs(outdir, calls, freqs, ubiq, resol, cross,
                          consolidated, diags, prefix=args.prefix)
    (outdir / "diagnostico_em.json").write_text(json.dumps(diags, indent=2, default=str))

    print("\nResumo das decisoes:")
    print(calls["decision"].value_counts().to_string())
    print("\nArquivos gerados:")
    for f in files:
        print("  ", f)
    return 0


def cmd_benchmark(args) -> int:
    """Valida o pipeline contra uma tabela de verdade conhecida."""
    from .validate import carregar_verdade, comparar_com_verdade, masking_experiment

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[1/4] lendo saidas das ferramentas...")
    raw = load_sources(args.inputs, tool_map=_parse_sheet_tools(args.sheet_tool),
                       sample_regex=args.sample_regex)
    if raw.empty:
        print("ERRO: nenhuma chamada reconhecida nos caminhos informados.", file=sys.stderr)
        return 1
    print(f"      {len(raw)} registros | {raw['sample'].nunique()} amostras "
          f"| ferramentas: {sorted(raw['tool'].dropna().unique())}")

    print("[2/4] rodando o pipeline...")
    cross = cross_tool(consolidate(raw))
    calls, freqs, diags = resolve_cohort(cross, min_samples_for_em=args.min_samples_em)

    print("[3/4] comparando com a verdade...")
    verdade = carregar_verdade(args.truth, sample_regex=args.sample_regex)
    print(f"      verdade disponivel para {len(verdade)} pares amostra x locus")
    det, por_locus = comparar_com_verdade(calls, cross, verdade)
    if det.empty:
        print("ERRO: nenhum par amostra x locus em comum. Confira --sample-regex.",
              file=sys.stderr)
        return 1

    resumo = por_locus.attrs["resumo"]
    mask, met = masking_experiment(cross, level=1)

    print("[4/4] gravando...")
    det.to_csv(outdir / "benchmark_detalhado.tsv", sep="\t", index=False)
    por_locus.to_csv(outdir / "benchmark_por_locus.tsv", sep="\t", index=False)
    calls.drop(columns=["candidates"], errors="ignore").to_csv(
        outdir / "benchmark_chamadas.tsv", sep="\t", index=False)
    (outdir / "benchmark_resumo.json").write_text(
        json.dumps({"verdade": resumo, "mascaramento": met}, indent=2, default=str))

    print("\nRESULTADO")
    print(f"  pares com verdade      : {resumo['n_com_verdade']}")
    print(f"  cobertura (deu chamada): {resumo['cobertura']:.1%}")
    print(f"  acuracia               : {resumo['acuracia']:.1%}")
    print(f"  verdade entre candidatos: {resumo['sensibilidade_candidatos']:.1%}")
    print(f"  verdade descartada      : {resumo['verdade_descartada']}  <- erro grave se > 0")
    print("\nPor locus:")
    print(por_locus.to_string(index=False))
    return 0


def _parse_sheet_tools(pairs) -> dict[str, str]:
    out = {}
    for p in pairs or []:
        if "=" not in p:
            raise SystemExit(f"--sheet-tool espera ABA=FERRAMENTA, recebi {p!r}")
        k, v = p.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="kirsolve",
        description="Resolucao de alelos KIR ambiguos (PING + kir-mapper)",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("input", nargs="+",
                        help="planilha .xlsx, arquivos csv/tsv ou diretorios com as saidas")
    common.add_argument("--sheet-tool", action="append", metavar="ABA=FERRAMENTA",
                        help="associa uma aba a uma ferramenta, ex.: Sheet5=kir-mapper")
    common.add_argument("--sample-regex", default=DEFAULT_SAMPLE_REGEX,
                        help="regex para extrair o ID curto da amostra")

    pi = sub.add_parser("inspect", parents=[common],
                        help="mostra como o arquivo foi interpretado")
    pi.set_defaults(func=cmd_inspect)

    pr = sub.add_parser("run", parents=[common], help="roda o pipeline completo")
    pr.add_argument("-o", "--outdir", default="resultados_kirsolve")
    pr.add_argument("--prefix", default="kirsolve")
    pr.add_argument("--priors", help="TSV/CSV com locus, allele, frequency")
    pr.add_argument("--use-ipd", action="store_true",
                    help="baixa o catalogo IPD-KIR e usa como prior/validacao")
    pr.add_argument("--prior-weight", type=float, default=1.0,
                    help="peso do prior de Dirichlet (equivale a N amostras)")
    pr.add_argument("--min-samples-em", type=int, default=20,
                    help="minimo de amostras informativas para confiar no EM")
    pr.set_defaults(func=cmd_run)

    pb = sub.add_parser("benchmark",
                        help="valida o pipeline contra uma tabela de verdade")
    pb.add_argument("inputs", nargs="+",
                    help="arquivos ou diretorios com saidas de PING e kir-mapper")
    pb.add_argument("--truth", required=True,
                    help="tabela de verdade (HPRC/Immuannot, IHWG, Sanger...)")
    pb.add_argument("-o", "--outdir", default="benchmark")
    pb.add_argument("--sheet-tool", action="append", metavar="ARQUIVO=FERRAMENTA")
    pb.add_argument("--sample-regex", default=DEFAULT_SAMPLE_REGEX)
    pb.add_argument("--min-samples-em", type=int, default=20)
    pb.set_defaults(func=cmd_benchmark)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
