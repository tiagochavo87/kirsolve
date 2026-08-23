"""Testes unitarios: rode com `pytest -q` a partir da raiz do projeto."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kirsolve.em import posterior_genotypes, run_em
from kirsolve.nomenclature import (
    Allele,
    canonical_locus,
    find_locus_in_text,
    make_genotype,
    merge_genotypes,
)
from kirsolve.parsing import parse_call_string, split_composite_locus
from kirsolve.resolve import ambiguity_profile, safe_resolution


# ------------------------------------------------------------ nomenclatura

def test_parse_alelo_completo():
    a = Allele.parse("KIR2DL1*0010101")
    assert a.locus == "KIR2DL1" and a.fields == ("001", "01", "01") and a.n_fields == 3


def test_parse_sem_prefixo_kir():
    assert str(Allele.parse("2DL3*008N")) == "KIR2DL3*008N"


def test_parse_sufixo_e_truncagem():
    a = Allele.parse("KIR3DL1*0150101")
    assert str(a.truncate(1)) == "KIR3DL1*015"
    assert str(a.truncate(2)) == "KIR3DL1*01501"


def test_tokens_especiais():
    assert Allele.parse("KIR2DS2*null").is_null
    assert not Allele.parse("KIR2DL4*unresolved").is_typed


def test_variantes_novas_do_ping():
    a = Allele.parse("KIR2DL3*0010101$E4_13.G^E4_15.A")
    assert a.is_novel_candidate and a.novel_variants == ("E4_13.G", "E4_15.A")
    assert str(a.bare()) == "KIR2DL3*0010101"


def test_2dl5_colapsa_para_locus_unico():
    assert canonical_locus("KIR2DL5A") == "KIR2DL5"
    assert Allele.parse("KIR2DL5B*00201").locus == "KIR2DL5"


def test_locus_em_texto_livre():
    assert find_locus_in_text("copy_number_3DL3") == "KIR3DL3"
    assert find_locus_in_text("2DL23 PING e 2DL2 -Mapper alelo 1") == "KIR2DL23"


def test_compatibilidade_por_prefixo():
    a, b = Allele.parse("KIR2DL1*001"), Allele.parse("KIR2DL1*0010101")
    assert a.compatible_with(b) and str(a.merge(b)) == "KIR2DL1*0010101"
    assert not a.compatible_with(Allele.parse("KIR2DL1*002"))


# ------------------------------------------------------------ parsing

def test_kir_mapper_ponto_e_virgula():
    txt = "KIR2DL3*0010101+KIR2DL3*008N;KIR2DL3*0010101+KIR2DL3*034"
    gts, _ = parse_call_string(txt, "KIR2DL3", copy_number=2)
    assert len(gts) == 2 and all(len(g) == 2 for g in gts)


def test_ping_separado_por_espaco():
    txt = "KIR3DL3*00901+KIR3DL3*00901 KIR3DL3*00901+KIR3DL3*00906"
    gts, _ = parse_call_string(txt, "KIR3DL3", copy_number=2)
    assert len(gts) == 2


def test_ambiguidade_com_barra_vira_produto():
    txt = "2DL3*0010101 / 2DL3*0010103 + 2DL3*008N / 2DL3*034"
    gts, _ = parse_call_string(txt, "KIR2DL3", copy_number=2)
    assert len(gts) == 4


def test_lista_de_alelos_sem_pareamento():
    txt = "KIR3DL3*00901 KIR3DL3*00906 KIR3DL3*081"
    gts, warns = parse_call_string(txt, "KIR3DL3", copy_number=2)
    assert len(gts) == 6  # combinacoes com repeticao de 3 alelos em 2 copias
    assert warns


def test_copy_number_zero_vira_null():
    gts, _ = parse_call_string("KIR2DS2*null+KIR2DS2*null", "KIR2DS2", copy_number=0)
    assert len(gts) == 1 and next(iter(gts))[0].is_null


def test_failed():
    gts, _ = parse_call_string("failed", "KIR2DL1")
    assert next(iter(gts))[0].digits == "FAILED"


def test_split_locus_composto():
    gts, _ = parse_call_string("KIR2DL2*00101+KIR2DL3*00101", "KIR2DL23")
    out = split_composite_locus(gts, "KIR2DL23", {"KIR2DL2": 1, "KIR2DL3": 1})
    assert str(next(iter(out["KIR2DL2"]))[0]) == "KIR2DL2*00101"
    assert str(next(iter(out["KIR2DL3"]))[0]) == "KIR2DL3*00101"


# ------------------------------------------------------------ resolucao

def test_merge_genotipos_mantem_maior_resolucao():
    g = make_genotype([Allele.parse("KIR2DL1*003"), Allele.parse("KIR2DL1*006")])
    h = make_genotype([Allele.parse("KIR2DL1*00302"), Allele.parse("KIR2DL1*00602")])
    m = merge_genotypes(g, h)
    assert {str(a) for a in m} == {"KIR2DL1*00302", "KIR2DL1*00602"}


def test_merge_incompativel_retorna_none():
    g = make_genotype([Allele.parse("KIR2DL1*003"), Allele.parse("KIR2DL1*006")])
    h = make_genotype([Allele.parse("KIR2DL1*001"), Allele.parse("KIR2DL1*002")])
    assert merge_genotypes(g, h) is None


def test_resolucao_segura_em_campo_1():
    cands = {
        make_genotype([Allele.parse("KIR2DP1*00201"), Allele.parse("KIR2DP1*017")]),
        make_genotype([Allele.parse("KIR2DP1*00204"), Allele.parse("KIR2DP1*017")]),
    }
    n, gt = safe_resolution(cands)
    assert n == 1 and {str(a) for a in gt} == {"KIR2DP1*002", "KIR2DP1*017"}
    assert ambiguity_profile(cands) == {1: 1, 2: 2, 3: 2}


def test_em_recupera_frequencias():
    A = Allele.parse("KIR2DL1*003")
    B = Allele.parse("KIR2DL1*006")
    homoA = make_genotype([A, A])
    het = make_genotype([A, B])
    # 80 amostras nao ambiguas: 60 AA, 20 AB -> p(A)=0.875
    sets = [{homoA}] * 60 + [{het}] * 20
    freqs, diag = run_em(sets, pseudocount=0.0)
    assert diag["converged"]
    assert abs(freqs[A] - 0.875) < 0.01


def test_em_desambigua_com_a_coorte():
    A = Allele.parse("KIR3DL1*001")
    B = Allele.parse("KIR3DL1*002")
    C = Allele.parse("KIR3DL1*099")  # raro
    certo = make_genotype([A, B])
    # 50 amostras resolvidas com A/B e 1 ambigua entre A+B e A+C
    sets = [{certo}] * 50 + [{certo, make_genotype([A, C])}]
    freqs, _ = run_em(sets, pseudocount=0.01)
    post = posterior_genotypes(sets[-1], freqs)
    assert post[0][0] == certo and post[0][1] > 0.95
