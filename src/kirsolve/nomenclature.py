"""
Nomenclatura KIR (IPD-KIR): parsing, normalizacao, niveis de resolucao
e harmonizacao entre os loci usados por PING e kir-mapper.

Formato oficial: KIR2DL1*0010101
  digitos 1-3  -> campo 1 (proteina distinta)
  digitos 4-5  -> campo 2 (substituicao sinonima)
  digitos 6-7  -> campo 3 (variacao nao-codificante)
  sufixo N/L/S/C/Q -> expressao

Tokens especiais encontrados nas saidas reais:
  *null        -> gene/haplotipo ausente (CN contribui com 0)
  *unresolved  -> ferramenta nao conseguiu discriminar nenhum alelo
  failed       -> falha tecnica na chamada
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import total_ordering

# ---------------------------------------------------------------- genes e loci

KIR_GENES = [
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
    "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
    "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DS1",
    "KIR2DP1", "KIR3DP1",
]

#: loci compostos que as ferramentas emitem quando nao separam os parentes proximos
COMPOSITE_LOCI = {
    "KIR2DL23":  ("KIR2DL2", "KIR2DL3"),
    "KIR2DS35":  ("KIR2DS3", "KIR2DS5"),
    "KIR3DL1S1": ("KIR3DL1", "KIR3DS1"),
    "KIR2DL5AB": ("KIR2DL5A", "KIR2DL5B"),
    "KIR2DL5":   ("KIR2DL5A", "KIR2DL5B"),
}

ALL_LOCI = set(KIR_GENES) | set(COMPOSITE_LOCI)

#: locus de comparacao usado no pipeline. 2DL5A/2DL5B sao colapsados porque
#: nem PING nem kir-mapper os separam de forma confiavel em leitura curta.
LOCUS_CANONICAL = {g: g for g in KIR_GENES}
LOCUS_CANONICAL.update({
    "KIR2DL5A": "KIR2DL5",
    "KIR2DL5B": "KIR2DL5",
    "KIR2DL5AB": "KIR2DL5",
    "KIR2DL5": "KIR2DL5",
})

#: loci finais reportados pelo pipeline
CANONICAL_LOCI = [
    "KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5", "KIR2DS3",
    "KIR2DS5", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1",
    "KIR3DS1", "KIR2DS1", "KIR2DS4", "KIR3DL2",
]

#: genes framework: esperados em praticamente todo haplotipo (CN tipico = 2)
FRAMEWORK_GENES = {"KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"}

EXPRESSION_SUFFIXES = "NLSCQ"

# ---------------------------------------------------------------- regex

# a ordem da alternancia importa: 3DL1S1 e 2DL5AB antes de 2DL23 antes de 2DL2
_GENE_CORE = r"(?:KIR)?([23]D[LSP](?:1S1|\d[AB]{1,2}|\d\d|\d))"
_ALLELE_RE = re.compile(
    rf"^\s*{_GENE_CORE}\s*\*\s*(\d{{3,9}})([{EXPRESSION_SUFFIXES}]?)\s*$", re.IGNORECASE
)
_SPECIAL_RE = re.compile(
    rf"^\s*{_GENE_CORE}\s*\*\s*(null|unresolved|unresolvable|failed|NA)\s*$", re.IGNORECASE
)
_GENE_ONLY_RE = re.compile(rf"^\s*{_GENE_CORE}\s*$", re.IGNORECASE)

NULL = "NULL"
UNRESOLVED = "UNRESOLVED"
FAILED = "FAILED"
SPECIAL_DIGITS = {NULL, UNRESOLVED, FAILED}

_SPECIAL_MAP = {
    "null": NULL,
    "unresolved": UNRESOLVED,
    "unresolvable": UNRESOLVED,
    "failed": FAILED,
    "na": FAILED,
}

#: celulas vazias / sem chamada
NULL_TOKENS = {"", "na", "n/a", "nan", "none", "null", "-", "--", ".", "?"}


class AlleleParseError(ValueError):
    pass


def _canon_gene_token(token: str) -> str | None:
    gene = "KIR" + token.upper()
    if gene in ALL_LOCI:
        return gene
    return None


def normalize_locus(text) -> str | None:
    """'2dl23' / 'KIR2DL5AB' -> nome de locus reconhecido, ou None."""
    if text is None:
        return None
    m = _GENE_ONLY_RE.match(str(text).strip())
    return _canon_gene_token(m.group(1)) if m else None


def find_locus_in_text(text) -> str | None:
    """Extrai o primeiro locus KIR de um texto livre, ex.: 'KIR2DL1_Copy_number'."""
    if text is None:
        return None
    best = None
    for m in re.finditer(_GENE_CORE, str(text), re.IGNORECASE):
        gene = _canon_gene_token(m.group(1))
        if gene and (best is None or len(gene) > len(best)):
            best = gene
    return best


def canonical_locus(locus: str) -> str:
    return LOCUS_CANONICAL.get(locus, locus)


def is_blank(text) -> bool:
    return text is None or str(text).strip().lower() in NULL_TOKENS


# ---------------------------------------------------------------- Allele


@total_ordering
@dataclass(frozen=True)
class Allele:
    """Alelo KIR normalizado, ou sentinela (NULL / UNRESOLVED / FAILED)."""

    locus: str
    digits: str
    suffix: str = ""
    #: variantes extras reportadas pelo PING apos '$' (candidato a alelo novo)
    novel_variants: tuple[str, ...] = ()

    # ------------------------------------------------------------ construcao
    @classmethod
    def parse(cls, text: str, locus_hint: str | None = None) -> "Allele":
        raw = str(text).strip()

        # anotacao de variantes novas do PING: KIR2DL3*0010101$E4_13.G^E4_15.A
        novel: tuple[str, ...] = ()
        if "$" in raw:
            raw, _, ann = raw.partition("$")
            novel = tuple(v for v in re.split(r"[\^]", ann) if v)
            raw = raw.strip()

        if raw.lower() in _SPECIAL_MAP:
            if not locus_hint:
                raise AlleleParseError(f"token {text!r} sem locus de contexto")
            return cls(canonical_locus(normalize_locus(locus_hint) or locus_hint),
                       _SPECIAL_MAP[raw.lower()])

        m = _SPECIAL_RE.match(raw)
        if m:
            locus = _canon_gene_token(m.group(1))
            if locus is None:
                raise AlleleParseError(f"locus desconhecido em {text!r}")
            return cls(canonical_locus(locus), _SPECIAL_MAP[m.group(2).lower()])

        m = _ALLELE_RE.match(raw)
        if m:
            locus = _canon_gene_token(m.group(1))
            digits, suffix = m.group(2), m.group(3).upper()
        else:
            m2 = re.fullmatch(rf"\s*(\d{{3,9}})([{EXPRESSION_SUFFIXES}]?)\s*", raw, re.IGNORECASE)
            if not (m2 and locus_hint):
                raise AlleleParseError(f"nao consegui interpretar o alelo: {text!r}")
            locus = normalize_locus(locus_hint) or locus_hint
            digits, suffix = m2.group(1), m2.group(2).upper()

        if locus is None:
            raise AlleleParseError(f"locus desconhecido em {text!r}")
        if len(digits) % 2 == 0:
            raise AlleleParseError(
                f"digitos invalidos ({len(digits)}) em {text!r}; esperado 3, 5, 7 ou 9"
            )
        return cls(canonical_locus(locus), digits, suffix, novel)

    @classmethod
    def try_parse(cls, text, locus_hint=None) -> "Allele | None":
        try:
            return cls.parse(text, locus_hint=locus_hint)
        except AlleleParseError:
            return None

    @classmethod
    def null(cls, locus: str) -> "Allele":
        return cls(canonical_locus(locus), NULL)

    @classmethod
    def unresolved(cls, locus: str) -> "Allele":
        return cls(canonical_locus(locus), UNRESOLVED)

    # ------------------------------------------------------------ propriedades
    @property
    def is_special(self) -> bool:
        return self.digits in SPECIAL_DIGITS

    @property
    def is_null(self) -> bool:
        return self.digits == NULL

    @property
    def is_typed(self) -> bool:
        return not self.is_special

    @property
    def is_novel_candidate(self) -> bool:
        return bool(self.novel_variants)

    @property
    def n_fields(self) -> int:
        return 0 if self.is_special else 1 + (len(self.digits) - 3) // 2

    @property
    def fields(self) -> tuple[str, ...]:
        if self.is_special:
            return ()
        d = self.digits
        return (d[:3],) + tuple(d[i:i + 2] for i in range(3, len(d), 2))

    # ------------------------------------------------------------ operacoes
    def truncate(self, n_fields: int) -> "Allele":
        """Reduz resolucao: truncate(1) -> KIR2DL1*001 (nivel proteina)."""
        if self.is_special:
            return self
        n = max(1, min(n_fields, self.n_fields))
        return Allele(self.locus, self.digits[:3 + 2 * (n - 1)], self.suffix)

    def bare(self) -> "Allele":
        """Remove anotacao de variantes novas (para agrupar em frequencias)."""
        return Allele(self.locus, self.digits, self.suffix)

    def is_prefix_of(self, other: "Allele") -> bool:
        if self.is_special or other.is_special:
            return self.digits == other.digits and self.locus == other.locus
        return self.locus == other.locus and other.digits.startswith(self.digits)

    def compatible_with(self, other: "Allele") -> bool:
        """Mesmo alelo em resolucoes possivelmente diferentes."""
        return self.is_prefix_of(other) or other.is_prefix_of(self)

    def merge(self, other: "Allele") -> "Allele":
        """Combina duas observacoes compativeis mantendo a maior resolucao."""
        if not self.compatible_with(other):
            raise ValueError(f"{self} e {other} nao sao compativeis")
        keep = self if self.n_fields >= other.n_fields else other
        return Allele(keep.locus, keep.digits, keep.suffix or other.suffix,
                      keep.novel_variants or other.novel_variants)

    # ------------------------------------------------------------ texto
    def __str__(self) -> str:
        if self.is_special:
            return f"{self.locus}*{self.digits.lower()}"
        tag = "$novel" if self.novel_variants else ""
        return f"{self.locus}*{self.digits}{self.suffix}{tag}"

    def __repr__(self) -> str:  # pragma: no cover
        return f"Allele({str(self)!r})"

    def _key(self):
        return (self.locus, self.is_special, self.digits, self.suffix)

    def __lt__(self, other: "Allele") -> bool:
        return self._key() < other._key()


# ---------------------------------------------------------------- genotipos

Genotype = tuple  # tuple[Allele, ...] ordenado


def make_genotype(alleles) -> Genotype:
    return tuple(sorted(alleles))


def genotype_to_gl(g: Genotype) -> str:
    return "+".join(str(a) for a in g)


def genotypes_to_gl(gs) -> str:
    """Conjunto de genotipos alternativos em GL string (separador '|')."""
    return "|".join(genotype_to_gl(g) for g in sorted(gs))


def truncate_genotype(g: Genotype, n_fields: int) -> Genotype:
    return make_genotype(a.truncate(n_fields) for a in g)


def genotypes_compatible(g: Genotype, h: Genotype) -> bool:
    """Existe pareamento 1-a-1 entre alelos compativeis (resolucoes distintas)."""
    if len(g) != len(h):
        return False
    return _match(list(g), list(h)) is not None


def merge_genotypes(g: Genotype, h: Genotype) -> Genotype | None:
    """Funde dois genotipos compativeis mantendo a maior resolucao de cada alelo."""
    pairing = _match(list(g), list(h))
    if pairing is None:
        return None
    return make_genotype(a.merge(b) for a, b in pairing)


def _match(g: list, h: list):
    """Pareamento perfeito por backtracking (n <= 4, custo irrelevante)."""
    if not g:
        return []
    a, rest = g[0], g[1:]
    for i, b in enumerate(h):
        if a.compatible_with(b):
            sub = _match(rest, h[:i] + h[i + 1:])
            if sub is not None:
                return [(a, b)] + sub
    return None
