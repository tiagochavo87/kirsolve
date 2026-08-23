"""
Priors de frequencia alelica para o EM.

Com coorte pequena (dezenas de amostras) o EM sozinho nao tem informacao para
desempatar ambiguidade. O prior de Dirichlet resolve isso trazendo frequencias
externas. Tres fontes, em ordem de preferencia:

  1. arquivo local TSV/CSV do usuario  (--priors freq.tsv)
     colunas: locus, allele, frequency
  2. lista de alelos do IPD-KIR (ANHIG/IPDKIR), usada como prior *uniforme
     informado*: so restringe o espaco aos alelos que existem no banco
  3. nenhum: prior uniforme com pseudocontagem

IMPORTANTE: use frequencias da populacao mais proxima da sua coorte. Frequencias
KIR variam muito entre grupos biogeograficos; um prior europeu aplicado a uma
coorte brasileira miscigenada enviesa a chamada.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import pandas as pd

from .nomenclature import Allele, canonical_locus

IPD_KIR_ALLELELIST_URL = (
    "https://raw.githubusercontent.com/ANHIG/IPDKIR/Latest/Allelelist.txt"
)


def load_prior_file(path: str | Path) -> dict[str, dict[Allele, float]]:
    """Le TSV/CSV com colunas locus, allele, frequency."""
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)
    cols = {c.lower().strip(): c for c in df.columns}
    need = {"allele", "frequency"}
    if not need <= set(cols):
        raise ValueError(f"{path} precisa das colunas 'allele' e 'frequency'")
    out: dict[str, dict[Allele, float]] = {}
    for _, r in df.iterrows():
        locus_hint = r[cols["locus"]] if "locus" in cols else None
        a = Allele.try_parse(r[cols["allele"]], locus_hint=locus_hint)
        if a is None or not a.is_typed:
            continue
        out.setdefault(a.locus, {})[a.bare()] = float(r[cols["frequency"]])
    # normaliza por locus
    for locus, d in out.items():
        tot = sum(d.values())
        if tot > 0:
            out[locus] = {a: v / tot for a, v in d.items()}
    return out


def parse_ipd_allelelist(text: str) -> dict[str, set[Allele]]:
    """Extrai o catalogo de alelos validos do Allelelist.txt do IPD-KIR."""
    catalog: dict[str, set[Allele]] = {}
    for line in io.StringIO(text):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for token in re.split(r"[,\t; ]+", line):
            a = Allele.try_parse(token)
            if a is not None and a.is_typed:
                catalog.setdefault(a.locus, set()).add(a.bare())
    return catalog


def fetch_ipd_catalog(cache: str | Path = "ipd_kir_allelelist.txt") -> dict[str, set[Allele]]:
    """
    Baixa (ou le do cache) a lista oficial de alelos do IPD-KIR.

    Requer rede na primeira execucao. Se falhar, retorna dicionario vazio e o
    pipeline segue sem restricao de catalogo.
    """
    cache = Path(cache)
    if cache.exists():
        return parse_ipd_allelelist(cache.read_text(encoding="utf-8", errors="ignore"))
    try:
        import urllib.request

        with urllib.request.urlopen(IPD_KIR_ALLELELIST_URL, timeout=60) as resp:
            text = resp.read().decode("utf-8", errors="ignore")
        cache.write_text(text, encoding="utf-8")
        return parse_ipd_allelelist(text)
    except Exception as exc:  # rede indisponivel, proxy, etc.
        print(f"[aviso] nao foi possivel obter o catalogo IPD-KIR ({exc}). "
              f"Seguindo sem restricao de catalogo.")
        return {}


def catalog_as_prior(
    catalog: dict[str, set[Allele]], observed: dict[str, set[Allele]] | None = None
) -> dict[str, dict[Allele, float]]:
    """Converte catalogo em prior uniforme sobre os alelos validos observados."""
    out = {}
    for locus, alleles in catalog.items():
        pool = alleles & observed[locus] if observed and locus in observed else alleles
        if pool:
            out[canonical_locus(locus)] = {a: 1.0 / len(pool) for a in pool}
    return out


def validate_against_catalog(
    calls_alleles: dict[str, set[Allele]], catalog: dict[str, set[Allele]]
) -> pd.DataFrame:
    """Lista alelos chamados que nao constam no catalogo IPD-KIR (possiveis novos)."""
    rows = []
    if not catalog:
        return pd.DataFrame(rows)
    for locus, alleles in calls_alleles.items():
        known = catalog.get(locus, set())
        if not known:
            continue
        for a in sorted(alleles):
            # compara por prefixo: 001 e valido se existe 0010101 no catalogo
            if not any(a.is_prefix_of(k) or k.is_prefix_of(a) for k in known):
                rows.append({"locus": locus, "allele": str(a),
                             "issue": "ausente do catalogo IPD-KIR"})
    return pd.DataFrame(rows)
