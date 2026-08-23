"""
Nucleo do pipeline: consolida evidencias, cruza PING x kir-mapper,
roda o EM na coorte e emite a chamada final com probabilidade posterior.

Ordem das operacoes (cada etapa so restringe o conjunto de candidatos):

  1. consolidate      - funde linhas duplicadas da mesma ferramenta (varias abas)
  2. cross_tool       - intersecao PING x kir-mapper por compatibilidade
  3. copy number      - ja aplicado no parsing
  4. safe_resolution  - nivel de campo em que a chamada e certa sem modelo
  5. EM + posterior   - desempate probabilistico entre os candidatos restantes
"""

from __future__ import annotations

from collections import Counter, defaultdict

import pandas as pd

from .em import posterior_genotypes, run_em
from .nomenclature import (
    Allele,
    CANONICAL_LOCI,
    Genotype,
    genotypes_to_gl,
    genotype_to_gl,
    merge_genotypes,
    truncate_genotype,
)

CONCORDANCE = {
    "concordant": "chamadas compativeis entre as ferramentas",
    "one_tool": "apenas uma ferramenta produziu chamada tipada",
    "discordant": "nenhum genotipo compativel entre as ferramentas",
    "no_call": "nenhuma ferramenta tipou o locus",
    "absent": "gene ausente (CN=0) nas duas ferramentas",
}


# --------------------------------------------------------------- consolidacao


def _merge_sets(a: set[Genotype], b: set[Genotype]) -> set[Genotype]:
    """Intersecao por compatibilidade, mantendo a maior resolucao de cada alelo."""
    if not a:
        return set(b)
    if not b:
        return set(a)
    out = set()
    for g in a:
        for h in b:
            m = merge_genotypes(g, h)
            if m is not None:
                out.add(m)
    return out


def _typed_only(gs: set[Genotype]) -> set[Genotype]:
    return {g for g in gs if any(x.is_typed for x in g)}


def consolidate(df: pd.DataFrame) -> pd.DataFrame:
    """Funde linhas da mesma (sample, tool, locus) vindas de abas diferentes."""
    rows = []
    for (sample, tool, locus), sub in df.groupby(["sample", "tool", "locus"], sort=False):
        typed_sets = [s for s in sub["genotypes"] if _typed_only(s)]
        acc: set[Genotype] = set()
        conflict = False
        for s in typed_sets:
            s = _typed_only(s)
            if not acc:
                acc = set(s)
                continue
            merged = _merge_sets(acc, s)
            if merged:
                acc = merged
            else:
                conflict = True
                acc |= s
        if not acc:
            # so sentinelas: mantem o mais informativo (null > unresolved)
            allsets = [s for s in sub["genotypes"] if s]
            acc = allsets[0] if allsets else set()

        cns = [c for c in sub["copy_number"].tolist() if pd.notna(c)]
        cn = int(Counter(cns).most_common(1)[0][0]) if cns else None
        cn_conflict = len(set(cns)) > 1

        rows.append({
            "sample": sample,
            "tool": tool,
            "locus": locus,
            "copy_number": cn,
            "cn_conflict": cn_conflict,
            "intra_tool_conflict": conflict,
            "genotypes": acc,
            "n_candidates": len(acc),
            "raw_calls": " || ".join(sorted({str(r) for r in sub["raw_call"] if str(r).strip()})),
            "sheets": ",".join(sorted(set(sub["source_sheet"]))),
            "warnings": "; ".join(sorted({w for w in sub["warnings"] if w})),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------- cruzamento


def cross_tool(df: pd.DataFrame, tools=("PING", "kir-mapper")) -> pd.DataFrame:
    """Cruza as ferramentas por (sample, locus) e classifica a concordancia."""
    rows = []
    for (sample, locus), sub in df.groupby(["sample", "locus"], sort=False):
        by_tool = {t: sub[sub["tool"] == t] for t in tools}
        sets, cns, extra = {}, {}, {}
        for t, s in by_tool.items():
            if len(s) == 0:
                sets[t], cns[t] = set(), None
                continue
            r = s.iloc[0]
            sets[t] = r["genotypes"]
            cns[t] = r["copy_number"]
            extra[t] = r

        typed = {t: _typed_only(g) for t, g in sets.items()}
        active = [t for t in tools if typed[t]]

        if len(active) == 2:
            inter = _merge_sets(typed[tools[0]], typed[tools[1]])
            if inter:
                status, candidates = "concordant", inter
            else:
                # Discordancia costuma ser artefato de resolucao, nao conflito real:
                # duas ferramentas podem apontar a mesma linhagem e divergir no 3o
                # campo. Tenta a intersecao em resolucao menor antes de desistir.
                status, candidates = None, None
                for lvl in (2, 1):
                    a = {truncate_genotype(g, lvl) for g in typed[tools[0]]}
                    b = {truncate_genotype(g, lvl) for g in typed[tools[1]]}
                    inter = _merge_sets(a, b)
                    if inter:
                        status, candidates = f"concordant_campo{lvl}", inter
                        break
                if candidates is None:
                    # conflito irreconciliavel: NUNCA unir (isso inflaria a
                    # ambiguidade acima da de cada ferramenta isolada).
                    # Mantem o conjunto menor e marca para inspecao manual.
                    smaller = min(tools, key=lambda t: len(typed[t]))
                    status, candidates = "discordant", typed[smaller]
        elif len(active) == 1:
            status, candidates = "one_tool", typed[active[0]]
        else:
            allcn = [c for c in cns.values() if c is not None]
            if allcn and all(c == 0 for c in allcn):
                status = "absent"
            else:
                status = "no_call"
            candidates = set()

        cn_vals = [c for c in cns.values() if c is not None]
        rows.append({
            "sample": sample,
            "locus": locus,
            "cn_PING": cns.get(tools[0]),
            "cn_kir_mapper": cns.get(tools[1]),
            "cn_agree": (len(set(cn_vals)) <= 1) if len(cn_vals) == 2 else None,
            "copy_number": int(max(set(cn_vals), key=cn_vals.count)) if cn_vals else None,
            "status": status,
            "n_candidates": len(candidates),
            "candidates": candidates,
            "n_PING": len(typed[tools[0]]),
            "n_kir_mapper": len(typed[tools[1]]),
            "raw_PING": extra.get(tools[0], {}).get("raw_calls", "") if tools[0] in extra else "",
            "raw_kir_mapper": extra.get(tools[1], {}).get("raw_calls", "") if tools[1] in extra else "",
        })
    out = pd.DataFrame(rows)
    return out.sort_values(["sample", "locus"]).reset_index(drop=True)


# --------------------------------------------------------------- resolucao


def safe_resolution(candidates: set[Genotype]) -> tuple[int, Genotype | None]:
    """
    Maior nivel de campo em que TODOS os candidatos colapsam num unico genotipo.

    Retorna (n_campos, genotipo) ou (0, None) se nem no nivel 1 ha unicidade.
    Este e o resultado 'certo', independente de qualquer modelo estatistico.
    """
    if not candidates:
        return 0, None
    for n in (3, 2, 1):
        collapsed = {truncate_genotype(g, n) for g in candidates}
        if len(collapsed) == 1:
            return n, next(iter(collapsed))
    return 0, None


def ambiguity_profile(candidates: set[Genotype]) -> dict[int, int]:
    """Numero de genotipos distintos em cada nivel de resolucao."""
    return {n: len({truncate_genotype(g, n) for g in candidates}) for n in (1, 2, 3)}


def resolve_cohort(
    cross: pd.DataFrame,
    priors: dict[str, dict[Allele, float]] | None = None,
    prior_weight: float = 1.0,
    top_k: int = 3,
    min_samples_for_em: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Roda o EM por locus e emite a chamada final.

    Retorna (chamadas, frequencias, diagnostico_por_locus).
    """
    priors = priors or {}
    freq_rows, call_rows, diags = [], [], {}
    freqs_by_locus: dict[str, dict[Allele, float]] = {}

    for locus in sorted(cross["locus"].unique(), key=_locus_order):
        sub = cross[cross["locus"] == locus]
        sets = [s for s in sub["candidates"] if s]
        freqs, diag = run_em(sets, prior=priors.get(locus), prior_weight=prior_weight)
        diag["em_reliable"] = diag["n_samples"] >= min_samples_for_em
        diags[locus] = diag
        freqs_by_locus[locus] = freqs
        for a, f in sorted(freqs.items(), key=lambda t: -t[1]):
            freq_rows.append({
                "locus": locus, "allele": str(a), "frequency": f,
                "n_samples_informative": diag["n_samples"],
                "em_reliable": diag["em_reliable"],
            })

    for _, r in cross.iterrows():
        cands = r["candidates"]
        n_safe, safe_gt = safe_resolution(cands)
        prof = ambiguity_profile(cands) if cands else {1: 0, 2: 0, 3: 0}
        post = posterior_genotypes(cands, freqs_by_locus.get(r["locus"], {}))
        best, best_p = (post[0] if post else (None, float("nan")))
        second_p = post[1][1] if len(post) > 1 else 0.0

        call_rows.append({
            "sample": r["sample"],
            "locus": r["locus"],
            "copy_number": r["copy_number"],
            "cn_agree": r["cn_agree"],
            "status": r["status"],
            "n_candidates": r["n_candidates"],
            "amb_field1": prof[1], "amb_field2": prof[2], "amb_field3": prof[3],
            "certain_level": n_safe,
            "certain_call": genotype_to_gl(safe_gt) if safe_gt else "",
            "best_call": genotype_to_gl(best) if best else "",
            "posterior": best_p,
            "margin": (best_p - second_p) if post else float("nan"),
            "em_reliable": diags.get(r["locus"], {}).get("em_reliable", False),
            "top_calls": " | ".join(
                f"{genotype_to_gl(g)}={p:.3f}" for g, p in post[:top_k]
            ),
            "gl_string": genotypes_to_gl(cands) if 0 < len(cands) <= 50 else "",
            "raw_PING": r["raw_PING"],
            "raw_kir_mapper": r["raw_kir_mapper"],
        })

    calls = pd.DataFrame(call_rows)
    calls["decision"] = calls.apply(_decision, axis=1)
    return calls, pd.DataFrame(freq_rows), diags


def _decision(row) -> str:
    if row["status"] in {"absent", "no_call"}:
        return row["status"]
    if row["status"] == "discordant":
        # As ferramentas nao se reconciliam nem em campo 1. Ficamos com o conjunto
        # menor para nao inflar a ambiguidade, mas isso NAO e resolucao: a escolha
        # foi por parcimonia, nao por evidencia. Exige inspecao do BAM.
        return "conflito_nao_resolvido"
    if row["n_candidates"] == 1:
        return "resolvido_unico"
    if row["certain_level"] == 3:
        return "resolvido_unico"
    if row["certain_level"] >= 1:
        return f"resolvido_campo{row['certain_level']}"
    if row["em_reliable"] and row["posterior"] >= 0.90 and row["margin"] >= 0.50:
        return "provavel_EM_alta"
    if row["em_reliable"] and row["posterior"] >= 0.70:
        return "provavel_EM_moderada"
    return "ambiguo"


def _locus_order(locus: str) -> int:
    try:
        return CANONICAL_LOCI.index(locus)
    except ValueError:
        return 999


# --------------------------------------------------------------- ubiquidade


def ubiquity_report(cross: pd.DataFrame) -> pd.DataFrame:
    """
    Identifica os 'alelos ubiquos': aqueles que aparecem em quase todo conjunto
    de candidatos mas quase nunca sao a chamada certa.

    - n_sets              : conjuntos ambiguos em que o alelo aparece
    - prop_sets           : proporcao dos conjuntos ambiguos do locus
    - n_certain           : conjuntos em que o alelo esta em TODOS os candidatos
    - ubiquity_index      : prop_sets * (1 - n_certain/n_sets); ~1 = puro ruido
                            de nomenclatura, ~0 = alelo realmente informativo
    """
    rows = []
    for locus, sub in cross.groupby("locus"):
        amb = [s for s in sub["candidates"] if len(s) > 1]
        if not amb:
            continue
        appear = Counter()
        certain = Counter()
        for s in amb:
            present = {a.bare() for g in s for a in g if a.is_typed}
            in_all = set.intersection(
                *[{a.bare() for a in g if a.is_typed} for g in s]
            ) if s else set()
            for a in present:
                appear[a] += 1
            for a in in_all:
                certain[a] += 1
        for a, n in appear.most_common():
            prop = n / len(amb)
            idx = prop * (1 - certain[a] / n)
            rows.append({
                "locus": locus,
                "allele": str(a),
                "n_ambiguous_sets": n,
                "total_ambiguous_sets": len(amb),
                "prop_sets": prop,
                "n_certain": certain[a],
                "ubiquity_index": idx,
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["locus", "ubiquity_index"], ascending=[True, False]).reset_index(drop=True)


def resolution_recommendation(cross: pd.DataFrame) -> pd.DataFrame:
    """
    Para cada locus, em qual nivel de campo vale a pena reportar:
    proporcao de amostras com chamada unica em 1, 2 e 3 campos.
    """
    rows = []
    for locus, sub in cross.groupby("locus"):
        sets = [s for s in sub["candidates"] if s]
        if not sets:
            continue
        n = len(sets)
        r = {"locus": locus, "n_samples": n}
        for lvl in (1, 2, 3):
            uniq = sum(1 for s in sets if len({truncate_genotype(g, lvl) for g in s}) == 1)
            r[f"prop_unique_field{lvl}"] = uniq / n
        r["recommended_level"] = max(
            (lvl for lvl in (3, 2, 1) if r[f"prop_unique_field{lvl}"] >= 0.80),
            default=0,
        )
        rows.append(r)
    return pd.DataFrame(rows).sort_values("locus").reset_index(drop=True)
