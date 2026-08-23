"""
Estimacao de frequencias alelicas por Expectation-Maximization (EM) a partir de
genotipos ambiguos, com copy number variavel.

Modelo
------
Para um locus com copy number c em uma amostra, o genotipo e um multiconjunto
de c alelos. Sob uniao aleatoria de gametas (equivalente a HWE generalizado
para c copias), a probabilidade de um multiconjunto g com contagens n_a e:

    P(g) = c! / prod_a(n_a!) * prod_a p_a^{n_a}

Cada amostra i contribui com um conjunto de candidatos C_i (a ambiguidade).
O EM alterna:

  E-step: w_{i,g} = P(g) / sum_{g' in C_i} P(g')
  M-step: p_a = (alpha_a + sum_i sum_g w_{i,g} n_a(g)) / (sum alpha + total copias)

`alpha` e um prior de Dirichlet, alimentado por frequencias externas
(IPD-KIR / Allele Frequency Net / literatura). Com coorte pequena o prior
domina, o que e o comportamento desejado: sem prior e sem N, o EM nao tem
informacao para desempatar.

Alelos sentinela (*null, *unresolved, *failed) nao entram no modelo.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np

from .nomenclature import Allele, Genotype


def genotype_log_prob(g: Genotype, logp: dict[Allele, float]) -> float:
    """log P(g) para o multiconjunto de alelos tipados de g."""
    typed = [a.bare() for a in g if a.is_typed]
    if not typed:
        return 0.0
    counts = Counter(typed)
    c = len(typed)
    lp = math.lgamma(c + 1)
    for a, n in counts.items():
        if a not in logp:
            return -math.inf
        lp += n * logp[a] - math.lgamma(n + 1)
    return lp


def _informative(candidates: set[Genotype]) -> bool:
    """Conjunto de candidatos e informativo se tem >=1 alelo tipado."""
    return any(a.is_typed for g in candidates for a in g)


def run_em(
    candidate_sets: list[set[Genotype]],
    prior: dict[Allele, float] | None = None,
    prior_weight: float = 1.0,
    pseudocount: float = 0.5,
    max_iter: int = 500,
    tol: float = 1e-9,
) -> tuple[dict[Allele, float], dict]:
    """
    Estima frequencias alelicas de um locus.

    Retorna (frequencias, diagnostico).
    """
    sets = [s for s in candidate_sets if s and _informative(s)]
    alleles = sorted({a.bare() for s in sets for g in s for a in g if a.is_typed})
    if not alleles:
        return {}, {"n_samples": 0, "n_alleles": 0, "iterations": 0,
                    "loglik": float("nan"), "converged": True}

    index = {a: i for i, a in enumerate(alleles)}
    k = len(alleles)

    # matriz esparsa: para cada amostra, lista de (vetor de contagens, coef multinomial)
    encoded = []
    for s in sets:
        rows = []
        for g in s:
            typed = [a.bare() for a in g if a.is_typed]
            if not typed:
                continue
            vec = np.zeros(k)
            for a in typed:
                vec[index[a]] += 1
            coef = math.lgamma(len(typed) + 1) - sum(
                math.lgamma(n + 1) for n in Counter(typed).values()
            )
            rows.append((vec, coef))
        if rows:
            encoded.append(rows)

    alpha = np.full(k, pseudocount)
    if prior:
        for a, f in prior.items():
            if a in index:
                alpha[index[a]] += prior_weight * max(f, 0.0)

    p = np.full(k, 1.0 / k)
    loglik_old = -np.inf
    converged = False
    it = 0

    for it in range(1, max_iter + 1):
        logp = np.log(np.clip(p, 1e-300, None))
        counts = np.zeros(k)
        loglik = 0.0
        for rows in encoded:
            lps = np.array([coef + vec @ logp for vec, coef in rows])
            m = lps.max()
            w = np.exp(lps - m)
            tot = w.sum()
            loglik += m + math.log(tot)
            w /= tot
            for wi, (vec, _) in zip(w, rows):
                counts += wi * vec
        new_p = counts + alpha
        new_p /= new_p.sum()
        if abs(loglik - loglik_old) < tol:
            p = new_p
            converged = True
            break
        p, loglik_old = new_p, loglik

    freqs = {a: float(p[i]) for i, a in enumerate(alleles)}
    diag = {
        "n_samples": len(encoded),
        "n_alleles": k,
        "iterations": it,
        "loglik": float(loglik_old),
        "converged": converged,
    }
    return freqs, diag


def posterior_genotypes(
    candidates: set[Genotype], freqs: dict[Allele, float], floor: float = 1e-6
) -> list[tuple[Genotype, float]]:
    """Ordena os candidatos por probabilidade posterior (soma 1 dentro do conjunto)."""
    if not candidates:
        return []
    logp = {a: math.log(max(f, floor)) for a, f in freqs.items()}
    scored = []
    for g in candidates:
        lp = genotype_log_prob(g, logp)
        if lp == -math.inf:
            lp = math.log(floor) * max(1, sum(1 for a in g if a.is_typed))
        scored.append((g, lp))
    m = max(lp for _, lp in scored)
    exp = [(g, math.exp(lp - m)) for g, lp in scored]
    tot = sum(v for _, v in exp)
    return sorted(((g, v / tot) for g, v in exp), key=lambda t: (-t[1], str(t[0])))


def bootstrap_frequencies(
    candidate_sets: list[set[Genotype]],
    n_boot: int = 200,
    seed: int = 0,
    **em_kwargs,
) -> dict[Allele, tuple[float, float]]:
    """IC 95% percentil para as frequencias, por reamostragem de amostras."""
    rng = np.random.default_rng(seed)
    sets = [s for s in candidate_sets if s and _informative(s)]
    if len(sets) < 5:
        return {}
    draws: dict[Allele, list[float]] = defaultdict(list)
    for _ in range(n_boot):
        idx = rng.integers(0, len(sets), len(sets))
        f, _ = run_em([sets[i] for i in idx], **em_kwargs)
        for a, v in f.items():
            draws[a].append(v)
    out = {}
    for a, vals in draws.items():
        arr = np.array(vals)
        out[a] = (float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)))
    return out
