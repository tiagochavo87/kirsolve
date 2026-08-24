"""
Parsing das strings de chamada de alelos produzidas por PING e kir-mapper.

Formatos suportados (todos observados nos dados reais):

kir-mapper (coluna *_Calls)
    "KIR2DL3*0010101+KIR2DL3*008N;KIR2DL3*0010101+KIR2DL3*034;..."
    genotipos alternativos separados por ';', alelos por '+'

PING (finalAlleleCalls / saida de chamada de variantes)
    "KIR3DL3*00901+KIR3DL3*00901 KIR3DL3*00901+KIR3DL3*00906 ..."
    genotipos alternativos separados por espaco
    "KIR2DL4*unresolved+KIR2DL4*unresolved"
    "failed"

PING (interallelecall, com variantes novas)
    "KIR2DL3*0010101$E4_13.G^E4_15.A+KIR2DL3*00201"
    lista de alelos individuais separados por espaco (nao pareados)

Planilha "limpa" (uma coluna por alelo, ambiguidade com '/')
    "2DL3*0010101 / 2DL3*0010103 / 2DL3*008N"
"""

from __future__ import annotations

import itertools
import re

from .nomenclature import (
    Allele,
    AlleleParseError,
    COMPOSITE_LOCI,
    Genotype,
    canonical_locus,
    is_blank,
    make_genotype,
)

#: separadores de genotipos alternativos, em ordem de precedencia
GENOTYPE_SEPARATORS = [";", "|", "  ", " "]


def _split_genotypes(text: str) -> list[str]:
    """Divide a string em genotipos alternativos sem quebrar dentro de um '+'."""
    # espacos ao redor de '+' e '/' sao cosmeticos e nao separam genotipos
    text = re.sub(r"\s*([+/])\s*", r"\1", text.strip())
    for sep in (";", "|"):
        if sep in text:
            return [t.strip() for t in text.split(sep) if t.strip()]
    # espaco so separa genotipos quando ha '+' (PING). Caso contrario e lista de alelos.
    parts = [t for t in re.split(r"\s+", text) if t]
    return parts


def parse_call_string(
    text,
    locus: str,
    copy_number: int | None = None,
    max_genotypes: int = 20000,
) -> tuple[set[Genotype], list[str]]:
    """
    Converte uma celula de chamada no conjunto de genotipos candidatos.

    Retorna (genotipos, avisos). Genotipos sao tuplas ordenadas de Allele,
    todos com o mesmo tamanho (o copy number observado, quando informado).
    """
    warnings: list[str] = []
    if is_blank(text):
        return set(), warnings

    raw = str(text).strip()
    if raw.lower() in {"failed", "fail"}:
        return {make_genotype([Allele(canonical_locus(locus), "FAILED")])}, warnings

    chunks = _split_genotypes(raw)
    genotypes: set[Genotype] = set()
    loose_alleles: list[Allele] = []

    for chunk in chunks:
        if "+" in chunk:
            options = []
            ok = True
            for token in chunk.split("+"):
                alts = _parse_allele_alternatives(token, locus, warnings)
                if not alts:
                    ok = False
                    break
                options.append(alts)
            if not ok:
                continue
            for combo in itertools.product(*options):
                genotypes.add(make_genotype(combo))
                if len(genotypes) > max_genotypes:
                    warnings.append(
                        f"{locus}: mais de {max_genotypes} genotipos candidatos, truncado"
                    )
                    return genotypes, warnings
        else:
            loose_alleles.extend(_parse_allele_alternatives(chunk, locus, warnings))

    # Lista de alelos soltos (interallelecall, ou coluna 'alelo N' da planilha limpa):
    # combina em genotipos de tamanho = copy number.
    if loose_alleles and not genotypes:
        uniq = sorted(set(loose_alleles))
        n = copy_number if copy_number and copy_number > 0 else 1
        if len(uniq) == n:
            genotypes.add(make_genotype(uniq))
        else:
            for combo in itertools.combinations_with_replacement(uniq, n):
                genotypes.add(make_genotype(combo))
            if len(uniq) > 1:
                warnings.append(
                    f"{locus}: {len(uniq)} alelos sem pareamento explicito, "
                    f"expandido para CN={n}"
                )

    if copy_number is not None and copy_number >= 0:
        genotypes = _enforce_copy_number(genotypes, locus, copy_number, warnings)

    return genotypes, warnings


def _parse_allele_alternatives(token: str, locus: str, warnings: list[str]) -> list[Allele]:
    """Um token pode conter ambiguidade de alelo com '/': 'A*001 / A*002'."""
    out: list[Allele] = []
    for piece in re.split(r"\s*/\s*", token.strip()):
        if is_blank(piece):
            continue
        try:
            out.append(Allele.parse(piece, locus_hint=locus))
        except AlleleParseError as exc:
            warnings.append(f"{locus}: {exc}")
    return out


def _enforce_copy_number(
    genotypes: set[Genotype], locus: str, cn: int, warnings: list[str]
) -> set[Genotype]:
    """
    Ajusta cada genotipo ao copy number observado.

    - alelos *null nao contam para o CN (representam ausencia do gene)
    - se sobram alelos demais, o genotipo e descartado (inconsistente)
    - se faltam, completa com *null quando CN==0, senao mantem e sinaliza
    """
    if not genotypes:
        return genotypes
    fixed: set[Genotype] = set()
    for g in genotypes:
        typed = [a for a in g if a.is_typed or a.digits == "UNRESOLVED"]
        if cn == 0:
            fixed.add(make_genotype([Allele.null(locus)]))
            continue
        if len(typed) == cn:
            fixed.add(make_genotype(typed))
        elif len(typed) < cn:
            # ferramenta reportou menos alelos que o CN: completa com 'unresolved'
            pad = [Allele.unresolved(locus)] * (cn - len(typed))
            fixed.add(make_genotype(typed + pad))
        else:
            # mais alelos que o CN: gera todas as sub-combinacoes plausiveis
            for combo in itertools.combinations(typed, cn):
                fixed.add(make_genotype(combo))
    if not fixed:
        warnings.append(f"{locus}: nenhum genotipo compativel com CN={cn}")
    return fixed


def split_composite_locus(
    genotypes: set[Genotype], composite: str, copy_numbers: dict[str, int] | None = None
) -> dict[str, set[Genotype]]:
    """
    Separa chamadas de um locus composto (2DL23, 2DS35, 3DL1S1) nos genes membros,
    usando o proprio nome do alelo. 'KIR2DL2*00101+KIR2DL3*00101' vira
    {KIR2DL2: {(00101,)}, KIR2DL3: {(00101,)}}.

    Quando o genotipo e inteiramente 'unresolved'/'null', propaga o sentinela
    para todos os genes membros.
    """
    members = COMPOSITE_LOCI.get(composite, (composite,))
    members = tuple(canonical_locus(m) for m in members)
    out: dict[str, set[Genotype]] = {m: set() for m in members}

    for g in genotypes:
        buckets: dict[str, list[Allele]] = {m: [] for m in members}
        unassigned: list[Allele] = []
        for a in g:
            target = canonical_locus(a.locus)
            if target in buckets:
                buckets[target].append(a)
            else:
                unassigned.append(a)
        for m in members:
            alleles = buckets[m]
            if alleles:
                out[m].add(make_genotype(alleles))
            elif unassigned:
                # sentinela nao atribuivel: replica em todos os membros
                out[m].add(make_genotype(
                    [Allele(m, a.digits) for a in unassigned]
                ))
    for m in members:
        if copy_numbers and m in copy_numbers and out[m]:
            out[m] = _enforce_copy_number(out[m], m, copy_numbers[m], [])
        if not out[m]:
            out[m] = set()
    return out
