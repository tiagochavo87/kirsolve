"""Relatorios de saida: Excel multi-aba, TSVs e resumo de QC."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .nomenclature import CANONICAL_LOCI, genotype_to_gl

SHEET_NAMES = {
    "calls": "01_chamadas_finais",
    "ambiguous": "02_ainda_ambiguos",
    "discordant": "03_discordancias",
    "freqs": "04_frequencias_EM",
    "ubiquity": "05_alelos_ubiquos",
    "resolution": "06_resolucao_recomendada",
    "cn_qc": "07_qc_copy_number",
    "novel": "08_candidatos_novos",
    "raw": "09_consolidado_por_ferramenta",
    "diag": "10_diagnostico_EM",
}


def copy_number_qc(cross: pd.DataFrame) -> pd.DataFrame:
    """Discrepancias de copy number entre as ferramentas e violacoes de framework."""
    from .nomenclature import FRAMEWORK_GENES

    rows = []
    for _, r in cross.iterrows():
        issues = []
        if r["cn_agree"] is False:
            issues.append(f"CN divergente: PING={r['cn_PING']} vs kir-mapper={r['cn_kir_mapper']}")
        cn = r["copy_number"]
        if r["locus"] in FRAMEWORK_GENES and cn is not None and cn == 0:
            issues.append("gene framework com CN=0 (improvavel; revisar cobertura)")
        if cn is not None and cn > 4:
            issues.append(f"CN={cn} atipicamente alto")
        if issues:
            rows.append({
                "sample": r["sample"], "locus": r["locus"],
                "cn_PING": r["cn_PING"], "cn_kir_mapper": r["cn_kir_mapper"],
                "issues": "; ".join(issues),
            })
    return pd.DataFrame(rows)


def novel_candidates(cross: pd.DataFrame) -> pd.DataFrame:
    """Alelos chamados com variantes adicionais (anotacao '$' do PING)."""
    rows = []
    for _, r in cross.iterrows():
        seen = set()
        for g in r["candidates"]:
            for a in g:
                if a.is_novel_candidate and str(a) not in seen:
                    seen.add(str(a))
                    rows.append({
                        "sample": r["sample"], "locus": r["locus"],
                        "closest_allele": str(a.bare()),
                        "extra_variants": ";".join(a.novel_variants),
                        "n_extra_variants": len(a.novel_variants),
                    })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.drop_duplicates().sort_values(
        ["n_extra_variants", "sample", "locus"], ascending=[False, True, True]
    ).reset_index(drop=True)


def summarize(calls: pd.DataFrame) -> pd.DataFrame:
    """Resumo por locus: quanto o pipeline resolveu."""
    rows = []
    for locus, sub in calls.groupby("locus"):
        typed = sub[~sub["status"].isin(["absent", "no_call"])]
        n = len(typed)
        rows.append({
            "locus": locus,
            "n_amostras": len(sub),
            "n_tipaveis": n,
            "resolvido_unico": int((typed["decision"] == "resolvido_unico").sum()),
            "resolvido_parcial": int(typed["decision"].str.startswith("resolvido_campo").sum()),
            "provavel_EM": int(typed["decision"].str.startswith("provavel_EM").sum()),
            "ambiguo": int((typed["decision"] == "ambiguo").sum()),
            "prop_resolvido": (
                (typed["decision"].str.startswith("resolvido").sum() / n) if n else float("nan")
            ),
            "mediana_candidatos": typed["n_candidates"].median() if n else float("nan"),
        })
    out = pd.DataFrame(rows)
    out["_o"] = out["locus"].map(lambda l: CANONICAL_LOCI.index(l) if l in CANONICAL_LOCI else 99)
    return out.sort_values("_o").drop(columns="_o").reset_index(drop=True)


def gl_string_matrix(calls: pd.DataFrame, column: str = "best_call") -> pd.DataFrame:
    """Matriz amostra x locus com a chamada final, pronta para analise a jusante."""
    piv = calls.pivot_table(index="sample", columns="locus", values=column,
                            aggfunc="first")
    cols = [c for c in CANONICAL_LOCI if c in piv.columns]
    return piv[cols].reset_index()


def write_outputs(
    outdir: str | Path,
    calls: pd.DataFrame,
    freqs: pd.DataFrame,
    ubiquity: pd.DataFrame,
    resolution: pd.DataFrame,
    cross: pd.DataFrame,
    consolidated: pd.DataFrame,
    diags: dict,
    prefix: str = "kirsolve",
) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ambiguous = calls[calls["decision"] == "ambiguo"].copy()
    discordant = calls[calls["status"] == "discordant"].copy()
    cnqc = copy_number_qc(cross)
    novel = novel_candidates(cross)
    resumo = summarize(calls)
    matrix = gl_string_matrix(calls)
    diag_df = pd.DataFrame(
        [{"locus": k, **v} for k, v in diags.items()]
    ) if diags else pd.DataFrame()

    raw = consolidated.drop(columns=["genotypes"], errors="ignore").copy()

    xlsx = outdir / f"{prefix}_relatorio.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        resumo.to_excel(xw, sheet_name="00_resumo_por_locus", index=False)
        _drop_obj(calls).to_excel(xw, sheet_name=SHEET_NAMES["calls"], index=False)
        matrix.to_excel(xw, sheet_name="01b_matriz_chamadas", index=False)
        _drop_obj(ambiguous).to_excel(xw, sheet_name=SHEET_NAMES["ambiguous"], index=False)
        _drop_obj(discordant).to_excel(xw, sheet_name=SHEET_NAMES["discordant"], index=False)
        freqs.to_excel(xw, sheet_name=SHEET_NAMES["freqs"], index=False)
        ubiquity.to_excel(xw, sheet_name=SHEET_NAMES["ubiquity"], index=False)
        resolution.to_excel(xw, sheet_name=SHEET_NAMES["resolution"], index=False)
        cnqc.to_excel(xw, sheet_name=SHEET_NAMES["cn_qc"], index=False)
        novel.to_excel(xw, sheet_name=SHEET_NAMES["novel"], index=False)
        raw.to_excel(xw, sheet_name=SHEET_NAMES["raw"], index=False)
        diag_df.to_excel(xw, sheet_name=SHEET_NAMES["diag"], index=False)

    written = [xlsx]
    for name, df in [
        ("chamadas_finais", _drop_obj(calls)),
        ("matriz_chamadas", matrix),
        ("frequencias_EM", freqs),
        ("alelos_ubiquos", ubiquity),
        ("qc_copy_number", cnqc),
    ]:
        p = outdir / f"{prefix}_{name}.tsv"
        df.to_csv(p, sep="\t", index=False)
        written.append(p)
    return written


def _drop_obj(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=[c for c in ("candidates", "genotypes") if c in df.columns])
