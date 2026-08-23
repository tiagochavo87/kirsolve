"""
Validacao do pipeline.

Tres estrategias que usam apenas os proprios dados, sem material externo:

1. `masking_experiment`  - validacao cruzada leave-one-out. Pega as chamadas que
   sairam resolvidas, apaga artificialmente a resolucao (simulando uma
   ferramenta que so chegou ao campo 1), estima as frequencias com TODAS as
   OUTRAS amostras e verifica se o modelo recupera a resposta certa.
   E a medida direta da acuracia da etapa de desempate estatistico.

2. `hardy_weinberg` - para loci com CN=2, testa se as frequencias genotipicas
   observadas batem com as esperadas sob equilibrio. Desvio forte indica
   erro sistematico de chamada (alelo nulo nao detectado, cross-alignment)
   ou estrutura populacional.

3. `compare_frequencies` - confronta as frequencias estimadas com uma
   referencia externa (AFND, literatura). Divergencia grande em alelos comuns
   e sinal de vies.

O que este modulo NAO substitui: verdade experimental. Ver README, secao de
validacao externa.
"""

from __future__ import annotations

import math
from collections import Counter

import pandas as pd

from .em import posterior_genotypes, run_em
from .nomenclature import (
    Allele,
    Genotype,
    genotype_to_gl,
    genotypes_compatible,
    make_genotype,
    truncate_genotype,
)


# ------------------------------------------------------------ 1. mascaramento


def _plausible_candidates(
    truth: Genotype, pool: set[Allele], level: int
) -> set[Genotype]:
    """
    Reconstroi a ambiguidade que existiria se a ferramenta so tivesse chegado
    ao nivel `level`: todos os genotipos do pool que truncam para a mesma
    assinatura que a verdade.
    """
    alvo = truncate_genotype(truth, level)
    n = len(truth)
    cands = set()
    pool_l = sorted(pool)
    if n == 1:
        for a in pool_l:
            if truncate_genotype((a,), level) == alvo:
                cands.add(make_genotype([a]))
    elif n == 2:
        for i, a in enumerate(pool_l):
            for b in pool_l[i:]:
                if truncate_genotype((a, b), level) == alvo:
                    cands.add(make_genotype([a, b]))
    else:
        import itertools
        for combo in itertools.combinations_with_replacement(pool_l, n):
            if truncate_genotype(combo, level) == alvo:
                cands.add(make_genotype(combo))
    cands.add(truth)
    return cands


def masking_experiment(
    cross: pd.DataFrame,
    level: int = 1,
    prior_weight: float = 1.0,
    priors: dict | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Validacao cruzada leave-one-out sobre as chamadas ja resolvidas.

    Para cada amostra x locus com genotipo unico:
      - reconstroi a ambiguidade que haveria no nivel `level`
      - estima frequencias com todas as OUTRAS amostras do locus
      - verifica se o genotipo verdadeiro fica em 1o lugar no posterior

    Retorna (tabela por caso, metricas agregadas).
    """
    priors = priors or {}
    linhas = []

    for locus, sub in cross.groupby("locus"):
        sets = list(sub["candidates"])
        # a verdade e comparada sem anotacao de variante nova: '$novel' descreve
        # o alinhamento, nao a identidade do alelo
        resolvidos = [
            (i, make_genotype(a.bare() for a in next(iter(s))))
            for i, s in enumerate(sets)
            if len(s) == 1 and all(a.is_typed for a in next(iter(s)))
        ]
        if not resolvidos:
            continue
        pool = {a.bare() for s in sets for g in s for a in g if a.is_typed}
        if len(pool) < 2:
            continue

        idx = sub.index.tolist()
        for pos, verdade in resolvidos:
            cands = _plausible_candidates(verdade, pool, level)
            if len(cands) < 2:
                continue  # nao havia ambiguidade a recuperar

            treino = [s for j, s in enumerate(sets) if j != pos]
            freqs, diag = run_em(treino, prior=priors.get(locus),
                                 prior_weight=prior_weight)
            post = posterior_genotypes(cands, freqs)
            if not post:
                continue
            melhor, p_melhor = post[0]
            # acerto = compatibilidade, nao igualdade literal. '*001' e '*0010101'
            # sao o mesmo alelo em resolucoes diferentes; contar como erro seria
            # medir a resolucao do banco, nao a acuracia do modelo.
            p_verdade = sum(p for g, p in post if genotypes_compatible(g, verdade))
            rank = next(
                (k for k, (g, _) in enumerate(post, 1) if genotypes_compatible(g, verdade)),
                None,
            )

            linhas.append({
                "sample": cross.loc[idx[pos], "sample"],
                "locus": locus,
                "verdade": genotype_to_gl(verdade),
                "predito": genotype_to_gl(melhor),
                "acertou": genotypes_compatible(melhor, verdade),
                "n_candidatos_simulados": len(cands),
                "posterior_predito": p_melhor,
                "posterior_verdade": p_verdade,
                "rank_verdade": rank,
                "n_treino_informativo": diag["n_samples"],
            })

    df = pd.DataFrame(linhas)
    if df.empty:
        return df, {"n_casos": 0, "acuracia": float("nan"),
                    "acuracia_acaso": float("nan"), "aviso":
                    "nenhum caso com ambiguidade simulavel"}

    acaso = float((1.0 / df["n_candidatos_simulados"]).mean())
    met = {
        "n_casos": len(df),
        "acuracia": float(df["acertou"].mean()),
        "acuracia_acaso": acaso,
        "ganho_sobre_acaso": float(df["acertou"].mean()) - acaso,
        "top2": float((df["rank_verdade"] <= 2).mean()),
        "mediana_candidatos": float(df["n_candidatos_simulados"].median()),
        "mediana_treino": float(df["n_treino_informativo"].median()),
    }
    return df, met


# ------------------------------------------------------------ 2. Hardy-Weinberg


def hardy_weinberg(cross: pd.DataFrame, min_amostras: int = 30) -> pd.DataFrame:
    """
    Teste qui-quadrado de equilibrio de Hardy-Weinberg por locus, usando apenas
    amostras com CN=2 e genotipo unico. Agrupa alelos raros para manter as
    contagens esperadas viaveis.
    """
    linhas = []
    for locus, sub in cross.groupby("locus"):
        gts = [
            next(iter(s)) for _, r in sub.iterrows()
            for s in [r["candidates"]]
            if r["copy_number"] == 2 and len(s) == 1
            and all(a.is_typed for a in next(iter(s)))
        ]
        n = len(gts)
        if n < min_amostras:
            linhas.append({"locus": locus, "n_amostras": n, "chi2": None,
                           "gl": None, "p_valor": None,
                           "conclusao": f"n insuficiente (< {min_amostras})"})
            continue

        cont = Counter()
        for g in gts:
            cont[g] += 1
        alelos = Counter()
        for g, c in cont.items():
            for a in g:
                alelos[a.bare()] += c
        total = 2 * n
        p = {a: c / total for a, c in alelos.items()}

        chi2, gl = 0.0, 0
        for g, obs in cont.items():
            a, b = g[0].bare(), g[1].bare()
            esp = n * (p[a] ** 2 if a == b else 2 * p[a] * p[b])
            if esp < 1:
                continue
            chi2 += (obs - esp) ** 2 / esp
            gl += 1
        gl = max(gl - len(p), 1)
        pval = _chi2_sf(chi2, gl)
        linhas.append({
            "locus": locus, "n_amostras": n, "n_alelos": len(p),
            "chi2": round(chi2, 3), "gl": gl,
            "p_valor": round(pval, 4),
            "conclusao": "desvio significativo (p < 0,05)" if pval < 0.05
                         else "compativel com HWE",
        })
    return pd.DataFrame(linhas)


def _chi2_sf(x: float, k: int) -> float:
    """Cauda superior da qui-quadrado, sem scipy."""
    if x <= 0:
        return 1.0
    return _gammaincc(k / 2.0, x / 2.0)


def _gammaincc(a: float, x: float) -> float:
    if x < a + 1:
        return 1.0 - _gammainc_series(a, x)
    return _gammainc_cf(a, x)


def _gammainc_series(a, x, itmax=300, eps=1e-12):
    ap, s, d = a, 1.0 / a, 1.0 / a
    for _ in range(itmax):
        ap += 1
        d *= x / ap
        s += d
        if abs(d) < abs(s) * eps:
            break
    return s * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gammainc_cf(a, x, itmax=300, eps=1e-12):
    tiny = 1e-300
    b, c, d = x + 1 - a, 1 / tiny, 1 / (x + 1 - a)
    h = d
    for i in range(1, itmax):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < eps:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


# ------------------------------------------------- 3. comparacao de frequencias


def compare_frequencies(
    freqs: pd.DataFrame, referencia: dict[str, dict[Allele, float]]
) -> pd.DataFrame:
    """Compara frequencias estimadas com uma referencia externa."""
    linhas = []
    for _, r in freqs.iterrows():
        ref = referencia.get(r["locus"], {})
        a = Allele.try_parse(r["allele"])
        if a is None:
            continue
        f_ref = next((v for k, v in ref.items() if k.compatible_with(a)), None)
        if f_ref is None:
            continue
        linhas.append({
            "locus": r["locus"], "allele": r["allele"],
            "freq_estimada": r["frequency"], "freq_referencia": f_ref,
            "diferenca": r["frequency"] - f_ref,
            "razao": (r["frequency"] / f_ref) if f_ref > 0 else float("inf"),
        })
    df = pd.DataFrame(linhas)
    if df.empty:
        return df
    return df.sort_values("diferenca", key=abs, ascending=False).reset_index(drop=True)


# ------------------------------------------------- 4. comparacao com a verdade


def carregar_verdade(path, sample_regex: str = r"^([A-Za-z]+\d+)") -> dict:
    """
    Le a tabela de verdade (HPRC/Immuannot, IHWG, Sanger...).

    Formatos aceitos, detectados automaticamente:
      sample, locus, genotype           -> 'KIR2DL1*0030201+KIR2DL1*0030201'
      sample, locus, allele_1, allele_2 -> um alelo por coluna
      matriz amostra x locus            -> uma chamada por celula

    Devolve {(sample, locus): Genotype}.
    """
    import re
    from pathlib import Path

    from .loaders import normalize_sample
    from .nomenclature import canonical_locus, find_locus_in_text, is_blank

    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, dtype=str)
    else:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(path, sep=sep, dtype=str)

    cols = {str(c).strip().lower(): c for c in df.columns}
    verdade: dict = {}

    col_amostra = next((cols[k] for k in ("sample", "amostra", "id", "sample_id")
                        if k in cols), df.columns[0])
    col_locus = next((cols[k] for k in ("locus", "gene") if k in cols), None)

    if col_locus is not None:
        col_gt = next((cols[k] for k in ("genotype", "genotipo", "call", "alleles")
                       if k in cols), None)
        cols_alelo = [cols[k] for k in cols
                      if re.fullmatch(r"allele_?\d|alelo_?\d", k)]
        for _, r in df.iterrows():
            locus = canonical_locus(find_locus_in_text(r[col_locus]) or "")
            if not locus:
                continue
            if col_gt is not None and not is_blank(r[col_gt]):
                partes = re.split(r"[+/;]", str(r[col_gt]))
            else:
                partes = [r[c] for c in cols_alelo if not is_blank(r[c])]
            alelos = [Allele.try_parse(x, locus_hint=locus) for x in partes]
            alelos = [a for a in alelos if a is not None and a.is_typed]
            if alelos:
                chave = (normalize_sample(r[col_amostra], sample_regex), locus)
                verdade[chave] = make_genotype(a.bare() for a in alelos)
    else:
        for _, r in df.iterrows():
            amostra = normalize_sample(r[col_amostra], sample_regex)
            for c in df.columns[1:]:
                locus = canonical_locus(find_locus_in_text(c) or "")
                if not locus or is_blank(r[c]):
                    continue
                partes = re.split(r"[+/;]", str(r[c]))
                alelos = [Allele.try_parse(x, locus_hint=locus) for x in partes]
                alelos = [a for a in alelos if a is not None and a.is_typed]
                if alelos:
                    verdade[(amostra, locus)] = make_genotype(a.bare() for a in alelos)
    return verdade


def comparar_com_verdade(calls: pd.DataFrame, cross: pd.DataFrame,
                         verdade: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Confronta as chamadas do pipeline com a verdade conhecida.

    A comparacao respeita a resolucao: uma chamada reportada em campo 2 e
    julgada em campo 2. Cobrar 3 campos de quem so afirmou 2 mediria a
    resolucao do sequenciamento, nao a acuracia do metodo.

    Resultados possiveis por caso:
      correto              - a chamada bate com a verdade no nivel reportado
      incorreto            - nao bate
      verdade_descartada   - a verdade nem sequer estava entre os candidatos
                             (erro grave: alguma etapa a eliminou)
      verdade_entre_cands  - estava entre os candidatos, mas nao foi a escolhida
      sem_chamada          - ambiguo, conflito ou ausente
    """
    cands = {(r["sample"], r["locus"]): r["candidates"] for _, r in cross.iterrows()}
    linhas = []

    for _, r in calls.iterrows():
        chave = (r["sample"], r["locus"])
        if chave not in verdade:
            continue
        real = verdade[chave]
        conjunto = cands.get(chave, set())
        na_lista = any(genotypes_compatible(
            make_genotype(a.bare() for a in g), real) for g in conjunto)

        reportado = r["certain_call"] or r["best_call"]
        nivel = int(r["certain_level"]) if r["certain_level"] else 3

        if r["decision"] in {"absent", "no_call", "ambiguo", "conflito_nao_resolvido"}:
            resultado = "sem_chamada"
            acerto = None
        else:
            g = _parse_gl(reportado)
            if g is None:
                resultado, acerto = "sem_chamada", None
            else:
                alvo = truncate_genotype(real, nivel) if nivel else real
                acerto = genotypes_compatible(g, alvo)
                if acerto:
                    resultado = "correto"
                elif not na_lista:
                    resultado = "verdade_descartada"
                else:
                    resultado = "verdade_entre_cands"

        linhas.append({
            "sample": r["sample"], "locus": r["locus"],
            "verdade": genotype_to_gl(real),
            "reportado": reportado,
            "nivel_reportado": nivel,
            "decision": r["decision"],
            "resultado": resultado,
            "acerto": acerto,
            "verdade_nos_candidatos": na_lista,
            "n_candidatos": len(conjunto),
        })

    det = pd.DataFrame(linhas)
    if det.empty:
        return det, pd.DataFrame()

    chamadas = det[det["acerto"].notna()]
    resumo = {
        "n_com_verdade": len(det),
        "n_com_chamada": len(chamadas),
        "cobertura": len(chamadas) / len(det),
        "acuracia": float(chamadas["acerto"].mean()) if len(chamadas) else float("nan"),
        "sensibilidade_candidatos": float(det["verdade_nos_candidatos"].mean()),
        "verdade_descartada": int((det["resultado"] == "verdade_descartada").sum()),
    }
    por_locus = (det.groupby("locus")
                 .agg(n=("resultado", "size"),
                      acuracia=("acerto", "mean"),
                      verdade_nos_candidatos=("verdade_nos_candidatos", "mean"))
                 .reset_index())
    por_locus.attrs["resumo"] = resumo
    return det, por_locus


def _parse_gl(texto: str):
    """'KIR2DL1*003+KIR2DL1*003' -> Genotype, ignorando anotacao $novel."""
    if not texto or not str(texto).strip():
        return None
    alelos = []
    for parte in str(texto).split("+"):
        a = Allele.try_parse(parte.replace("$novel", "").strip())
        if a is None or not a.is_typed:
            return None
        alelos.append(a.bare())
    return make_genotype(alelos) if alelos else None
