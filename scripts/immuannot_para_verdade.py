#!/usr/bin/env python3
"""
Converte as anotacoes Immuannot (GTF por haplotipo, Zenodo 8372992) na tabela
de verdade que o `kirsolve benchmark` espera.

Entrada : diretorio com AMOSTRA.1.gtf.gz e AMOSTRA.2.gtf.gz
Saida   : TSV com sample, locus, genotype, e colunas de qualidade

Como funciona
-------------
Cada arquivo e UM haplotipo. O genotipo da amostra e a uniao dos dois:
o alelo do haplotipo 1 mais o do haplotipo 2. Como a montagem e fasada, a
verdade ja vem resolvida - e por isso que ela serve de padrao-ouro.

O campo `template_distance` diz quantas bases a montagem difere do alelo de
referencia mais proximo:

    0  -> o alelo esta no IPD-KIR exatamente como anotado
    >0 -> a montagem carrega um alelo NOVO, e `template_allele` e apenas o
          vizinho mais proximo

Isso e decisivo na hora de ler o benchmark. Se a verdade for um alelo novo,
nem o PING nem o kir-mapper tem como acertar: o alelo nao existe no banco
deles. Contar isso como erro do kirsolve seria injusto. Por padrao o script
marca esses casos numa coluna separada em vez de descarta-los, para voce
poder calcular a acuracia com e sem eles.

Uso:
    python3 immuannot_para_verdade.py /dados/verdade/hprc -o /dados/verdade.tsv
    python3 immuannot_para_verdade.py /dados/verdade/hprc -o v.tsv --apenas-exatos
    python3 immuannot_para_verdade.py /dados/verdade/hprc -o v.tsv --amostras /dados/amostras.txt
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

ATRIBUTO = re.compile(r'(\w+)\s+"?([^";]+)"?\s*;')
NOME_ARQUIVO = re.compile(r"^(?P<amostra>.+)\.(?P<hap>[12])\.gtf(\.gz)?$")

#: loci que o kirsolve reporta. 2DL5A e 2DL5B sao colapsados em 2DL5, porque
#: nem PING nem kir-mapper os separam de forma confiavel em leitura curta.
COLAPSAR = {"KIR2DL5A": "KIR2DL5", "KIR2DL5B": "KIR2DL5"}

#: limiar entre "mesmo alelo, banco diferente" e "alelo novo de verdade".
#: Calibrado na distribuicao real do HPRC, que e bimodal com vao entre 3 e 8.
LIMIAR_NOVO = 5


def classificar(distancia_max: int) -> str:
    if distancia_max == 0:
        return "exata"
    if distancia_max <= LIMIAR_NOVO:
        return "quase_exata"
    return "alelo_novo"


def ler_haplotipo(path: Path) -> list[dict]:
    """Extrai as linhas 'gene' de um GTF de haplotipo."""
    registros = []
    abrir = gzip.open if path.suffix == ".gz" else open
    with abrir(path, "rt", errors="ignore") as fh:
        for linha in fh:
            if linha.startswith("#"):
                continue
            campos = linha.rstrip("\n").split("\t")
            if len(campos) < 9 or campos[2].strip() != "gene":
                continue
            attrs = dict(ATRIBUTO.findall(campos[8]))
            gene = attrs.get("gene_name", "")
            if not gene.startswith("KIR"):
                continue
            alelo = attrs.get("template_allele", "")
            if not alelo:
                continue
            try:
                dist = int(attrs.get("template_distance", "0"))
            except ValueError:
                dist = 0
            registros.append({
                "locus": COLAPSAR.get(gene, gene),
                "gene_original": gene,
                "alelo": alelo.strip(),
                "distancia": dist,
            })
    return registros


def montar(diretorio: Path, amostras: set[str] | None) -> tuple[list[dict], list[str]]:
    por_amostra: dict[str, dict[int, list[dict]]] = defaultdict(dict)
    avisos: list[str] = []

    for arq in sorted(diretorio.glob("*.gtf.gz")) + sorted(diretorio.glob("*.gtf")):
        m = NOME_ARQUIVO.match(arq.name)
        if not m:
            avisos.append(f"nome fora do padrao, ignorado: {arq.name}")
            continue
        amostra = m.group("amostra")
        if amostras and amostra not in amostras:
            continue
        por_amostra[amostra][int(m.group("hap"))] = ler_haplotipo(arq)

    linhas = []
    for amostra in sorted(por_amostra):
        haps = por_amostra[amostra]
        if set(haps) != {1, 2}:
            avisos.append(
                f"{amostra}: haplotipos incompletos ({sorted(haps)}), amostra pulada"
            )
            continue

        loci = {r["locus"] for h in haps.values() for r in h}
        for locus in sorted(loci):
            alelos, distancias, genes = [], [], []
            for hap in (1, 2):
                for r in haps[hap]:
                    if r["locus"] == locus:
                        alelos.append(r["alelo"])
                        distancias.append(r["distancia"])
                        genes.append(r["gene_original"])

            novos = sum(1 for d in distancias if d > 0)
            linhas.append({
                "sample": amostra,
                "locus": locus,
                "genotype": "+".join(sorted(alelos)),
                "copy_number": len(alelos),
                "genes_originais": ",".join(sorted(set(genes))),
                "n_alelos_novos": novos,
                "distancia_max": max(distancias) if distancias else 0,
                "classe": classificar(max(distancias) if distancias else 0),
                "verdade_exata": novos == 0,
            })
    return linhas, avisos


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diretorio", help="pasta com os GTF do Immuannot")
    ap.add_argument("-o", "--saida", default="verdade_hprc.tsv")
    ap.add_argument("--amostras", help="arquivo com um ID por linha (filtra)")
    ap.add_argument("--apenas-exatos", action="store_true",
                    help="mantem so os loci com distancia zero")
    ap.add_argument("--sem-novos", action="store_true",
                    help="descarta so os alelos genuinamente novos (distancia > %d), "
                         "mantendo os quase-exatos" % LIMIAR_NOVO)
    args = ap.parse_args()

    filtro = None
    if args.amostras:
        filtro = {l.strip() for l in Path(args.amostras).read_text().splitlines()
                  if l.strip() and not l.startswith("#")}

    linhas, avisos = montar(Path(args.diretorio), filtro)
    for a in avisos:
        print(f"[aviso] {a}", file=sys.stderr)
    if not linhas:
        print("ERRO: nenhuma anotacao KIR encontrada.", file=sys.stderr)
        return 1

    if args.apenas_exatos:
        antes = len(linhas)
        linhas = [l for l in linhas if l["classe"] == "exata"]
        print(f"filtrados {antes - len(linhas)} loci (mantidos so os de distancia zero)")
    elif args.sem_novos:
        antes = len(linhas)
        linhas = [l for l in linhas if l["classe"] != "alelo_novo"]
        print(f"filtrados {antes - len(linhas)} loci com alelo novo")

    colunas = ["sample", "locus", "genotype", "copy_number", "genes_originais",
               "n_alelos_novos", "distancia_max", "classe", "verdade_exata"]
    with open(args.saida, "w") as fh:
        fh.write("\t".join(colunas) + "\n")
        for l in linhas:
            fh.write("\t".join(str(l[c]) for c in colunas) + "\n")

    amostras = {l["sample"] for l in linhas}
    from collections import Counter
    cont = Counter(l["classe"] for l in linhas)
    n = len(linhas)
    print(f"\n{args.saida}")
    print(f"  {len(amostras)} amostras | {n} pares amostra x locus\n")
    rotulos = {
        "exata":       "distancia 0     - alelo identico ao do IPD-KIR",
        "quase_exata": f"distancia 1-{LIMIAR_NOVO}   - versao de banco ou variacao nao codificante",
        "alelo_novo":  f"distancia >{LIMIAR_NOVO}    - alelo novo; as ferramentas nao podem acertar",
    }
    for classe in ("exata", "quase_exata", "alelo_novo"):
        c = cont.get(classe, 0)
        print(f"  {c:4d} ({c/n:5.1%})  {rotulos[classe]}")
    aproveitavel = cont.get("exata", 0) + cont.get("quase_exata", 0)
    print(f"\n  Base utilizavel para acuracia: {aproveitavel} loci ({aproveitavel/n:.1%}).")
    print("  Rode tambem com --sem-novos e compare os dois numeros.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
