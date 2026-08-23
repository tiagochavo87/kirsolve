#!/usr/bin/env python3
"""
Gera um conjunto de exemplo sintetico, para que o repositorio possa ser
executado e testado sem nenhum dado de paciente.

As chamadas imitam os formatos reais de PING e kir-mapper, incluindo os casos
que mais quebram parser: loci compostos, tokens `unresolved`/`null`/`failed`,
anotacao `$` de variantes novas, cabecalho deslocado e ambiguidade com `/`.

Os genotipos sao sorteados de um conjunto fixo de alelos com frequencias
arbitrarias. Nao representam nenhuma populacao real e nao devem ser usados
para inferencia — servem apenas para exercitar o codigo.

Uso:
    python3 scripts/gerar_exemplo.py exemplo/ --amostras 60 --semente 42
    python -m kirsolve run exemplo/planilha_exemplo.xlsx -o /tmp/saida
    python -m kirsolve benchmark exemplo/ping exemplo/kirmapper \\
        --truth exemplo/verdade.tsv -o /tmp/bench
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd

#: alelos por locus, com peso relativo. Valores inventados.
CATALOGO = {
    "KIR3DL3": [("00901", 5), ("00201", 3), ("0080101", 2), ("00906", 1)],
    "KIR2DL3": [("0010101", 6), ("00201", 3), ("006", 1)],
    "KIR2DL2": [("0010102", 4), ("00301", 2)],
    "KIR2DL1": [("00302", 5), ("0040109", 3), ("00303", 2), ("00201", 2)],
    "KIR2DL4": [("0080101", 4), ("0010201", 3), ("0050107", 2)],
    "KIR3DL1": [("0010103", 4), ("0020103", 3), ("01502", 2)],
    "KIR3DL2": [("0010101", 4), ("0020101", 3), ("00701", 2)],
    "KIR2DS4": [("0010101", 4), ("0030101", 3)],
    "KIR2DP1": [("0010201", 4), ("007", 2)],
    "KIR2DS2": [("0010124", 3), ("00101", 2)],
}

#: probabilidade de o gene estar ausente (CN=0) na amostra
P_AUSENTE = {"KIR2DL2": 0.45, "KIR2DS2": 0.45, "KIR2DP1": 0.15}

FRAMEWORK = {"KIR3DL3", "KIR2DL4", "KIR3DL2"}


def sortear(rng: random.Random, locus: str) -> list[str]:
    if locus not in FRAMEWORK and rng.random() < P_AUSENTE.get(locus, 0.05):
        return []
    alelos = [a for a, _ in CATALOGO[locus]]
    pesos = [w for _, w in CATALOGO[locus]]
    return [f"{locus}*{a}" for a in rng.choices(alelos, weights=pesos, k=2)]


def truncar(alelo: str, n: int) -> str:
    locus, digitos = alelo.split("*")
    manter = 3 + 2 * (n - 1)
    return f"{locus}*{digitos[:manter]}"


def gerar(destino: Path, n_amostras: int, semente: int) -> None:
    rng = random.Random(semente)
    destino.mkdir(parents=True, exist_ok=True)
    (destino / "ping").mkdir(exist_ok=True)
    (destino / "kirmapper").mkdir(exist_ok=True)

    ids = [f"SIM{i:03d}" for i in range(1, n_amostras + 1)]
    loci = list(CATALOGO)
    verdade, ping_calls, ping_cn, km_linhas = [], [], [], []

    for amostra in ids:
        linha_ping, linha_cn = {"sample": amostra}, {"sample": amostra}
        linha_km = {"Sample": amostra}

        for locus in loci:
            gt = sortear(rng, locus)
            cn = len(gt)
            linha_cn[locus] = cn
            linha_km[f"{locus}_Copy_number"] = cn

            if cn == 0:
                linha_ping[locus] = f"{locus}*null+{locus}*null"
                linha_km[f"{locus}_Calls"] = f"{locus}*null"
                continue

            verdade.append({"sample": amostra, "locus": locus,
                            "genotype": "+".join(sorted(gt))})

            # PING: resolucao menor + alternativas separadas por espaco.
            # Ocasionalmente falha ou anota variante nova.
            sorte = rng.random()
            if sorte < 0.06:
                linha_ping[locus] = "failed"
            elif sorte < 0.12:
                linha_ping[locus] = f"{locus}*unresolved+{locus}*unresolved"
            else:
                base = "+".join(truncar(a, 2) for a in sorted(gt))
                alts = [base]
                if rng.random() < 0.35:
                    outro = rng.choice([a for a, _ in CATALOGO[locus]])
                    alts.append("+".join(sorted([truncar(gt[0], 2),
                                                 truncar(f"{locus}*{outro}", 2)])))
                if rng.random() < 0.15:
                    alts[0] += "$E4_13.G^E4_15.A"
                linha_ping[locus] = " ".join(dict.fromkeys(alts))

            # kir-mapper: resolucao cheia, alternativas com ';'
            if rng.random() < 0.10:
                linha_km[f"{locus}_Calls"] = f"{locus}*unresolved"
            else:
                cands = ["+".join(sorted(gt))]
                if rng.random() < 0.30:
                    outro = rng.choice([a for a, _ in CATALOGO[locus]])
                    cands.append("+".join(sorted([gt[0], f"{locus}*{outro}"])))
                linha_km[f"{locus}_Calls"] = ";".join(dict.fromkeys(cands))

        ping_calls.append(linha_ping)
        ping_cn.append(linha_cn)
        km_linhas.append(linha_km)

    df_ping = pd.DataFrame(ping_calls)
    df_cn = pd.DataFrame(ping_cn)
    df_km = pd.DataFrame(km_linhas)

    # PING grava sem cabecalho na coluna do ID: reproduz o deslocamento real
    df_ping.rename(columns={"sample": ""}).to_csv(
        destino / "ping" / "finalAlleleCalls.csv", index=False)
    df_cn.rename(columns={"sample": ""}).to_csv(
        destino / "ping" / "manualCopyNumberFrame.csv", index=False)
    df_km.to_csv(destino / "kirmapper" / "kir-mapper_genotype_calls.tsv",
                 sep="\t", index=False)

    pd.DataFrame(verdade).to_csv(destino / "verdade.tsv", sep="\t", index=False)

    # planilha multi-aba, no formato que o pipeline recebeu originalmente
    with pd.ExcelWriter(destino / "planilha_exemplo.xlsx", engine="openpyxl") as xw:
        df_ping.to_excel(xw, sheet_name="saida PING", index=False)
        df_cn.to_excel(xw, sheet_name="copynumber PING", index=False)
        df_km.to_excel(xw, sheet_name="Sheet5", index=False)

    print(f"gerado em {destino}/")
    print(f"  {len(ids)} amostras x {len(loci)} loci")
    print(f"  verdade.tsv: {len(verdade)} pares amostra x locus")
    print(f"  ping/finalAlleleCalls.csv, ping/manualCopyNumberFrame.csv")
    print(f"  kirmapper/kir-mapper_genotype_calls.tsv")
    print(f"  planilha_exemplo.xlsx (3 abas)")
    print("\nDados sinteticos. Nao representam populacao real.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("destino", nargs="?", default="exemplo")
    ap.add_argument("--amostras", type=int, default=60)
    ap.add_argument("--semente", type=int, default=42)
    args = ap.parse_args()
    gerar(Path(args.destino), args.amostras, args.semente)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
