"""
Leitura das planilhas/tabelas de saida do PING e do kir-mapper.

O modulo detecta automaticamente qual dos layouts abaixo cada aba usa:

  A) "calls"  : Sample | software | <LOCUS>_Copy_number | <LOCUS>_Calls | ...
  B) "matrix" : Sample | <LOCUS> | <LOCUS> | ...   (uma chamada por celula)
  C) "cn"     : Sample | <LOCUS> | ...             (apenas inteiros = copy number)
  D) "wide"   : Sample | tool | copy_number_<X> | "<X> alelo 1" | "<X> alelo 2" | ...

O resultado e sempre a mesma tabela longa (`load_workbook` -> DataFrame):
    sample, tool, locus, copy_number, raw_call, n_candidates, genotypes, warnings
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .nomenclature import (
    COMPOSITE_LOCI,
    canonical_locus,
    find_locus_in_text,
    is_blank,
    normalize_locus,
)
from .parsing import parse_call_string, split_composite_locus

#: como reduzir IDs longos do PING (AMOSTRA001_23154FL-28Q2-01-140_S140_L004_R -> AMOSTRA001)
DEFAULT_SAMPLE_REGEX = r"^([A-Za-z]+\d+)"

TOOL_ALIASES = {
    "ping": "PING",
    "kir-mapper": "kir-mapper",
    "kirmapper": "kir-mapper",
    "kir_mapper": "kir-mapper",
    "mapper": "kir-mapper",
}


def normalize_tool(text) -> str | None:
    if is_blank(text):
        return None
    key = re.sub(r"[^a-z\-_]", "", str(text).strip().lower())
    return TOOL_ALIASES.get(key, str(text).strip())


def normalize_sample(text, pattern: str = DEFAULT_SAMPLE_REGEX) -> str:
    s = str(text).strip()
    m = re.match(pattern, s)
    return m.group(1) if m else s


@dataclass
class Record:
    sample: str
    tool: str
    locus: str
    copy_number: int | None
    raw_call: str
    genotypes: set = field(default_factory=set)
    warnings: list = field(default_factory=list)
    source_sheet: str = ""


# ------------------------------------------------------------------ deteccao


def _header_row(df: pd.DataFrame, max_scan: int = 5) -> int | None:
    """Encontra a linha que funciona como cabecalho (a que contem mais loci KIR)."""
    best, best_n = None, 0
    for i in range(min(max_scan, len(df))):
        n = sum(1 for v in df.iloc[i] if find_locus_in_text(v))
        if n > best_n:
            best, best_n = i, n
    return best if best_n >= 3 else None


def header_offset(df: pd.DataFrame, header_idx: int) -> int:
    """
    Algumas abas vem sem a celula vazia do canto superior esquerdo: o cabecalho
    comeca ja com um locus na coluna 0, enquanto a coluna 0 dos dados contem o
    ID da amostra. Nesse caso os rotulos estao deslocados +1 em relacao aos dados.
    """
    header = df.iloc[header_idx]
    if find_locus_in_text(header.iloc[0]) is None:
        return 0
    body = df.iloc[header_idx + 1:, 0].dropna()
    if body.empty:
        return 0
    # se a coluna 0 dos dados nao contem chamadas do locus anunciado, houve deslocamento
    looks_like_call = body.astype(str).str.contains(r"\*", regex=True).mean()
    return 0 if looks_like_call > 0.5 else 1


def detect_layout(df: pd.DataFrame, header_idx: int) -> str:
    header = [str(v) if pd.notna(v) else "" for v in df.iloc[header_idx]]
    joined = " ".join(header).lower()
    if "_calls" in joined or "calls" in joined and "copy_number" in joined:
        return "calls"
    if "alelo" in joined or "allele " in joined:
        return "wide"
    if "copy_number" in joined or "copy number" in joined:
        return "wide" if "alelo" in joined else "cn"
    body = df.iloc[header_idx + 1:, 1:]
    vals = [v for v in body.to_numpy().ravel() if pd.notna(v)]
    if vals and all(re.fullmatch(r"\s*\d+(\.0)?\s*", str(v)) for v in vals):
        return "cn"
    return "matrix"


# ------------------------------------------------------------------ leitores


def _read_calls_layout(df, header_idx, sheet, default_tool, sample_regex) -> list[Record]:
    header = [str(v) if pd.notna(v) else "" for v in df.iloc[header_idx]]
    cn_cols, call_cols = {}, {}
    tool_col = None
    for j, h in enumerate(header):
        hl = h.lower()
        if hl in {"software", "tool", "ping or kir-mapper", "programa"}:
            tool_col = j
            continue
        locus = find_locus_in_text(h)
        if not locus:
            continue
        if "copy" in hl:
            cn_cols[locus] = j
        elif "call" in hl or "alel" in hl or "allele" in hl:
            call_cols[locus] = j

    out: list[Record] = []
    for i in range(header_idx + 1, len(df)):
        row = df.iloc[i]
        if is_blank(row.iloc[0]):
            continue
        sample = normalize_sample(row.iloc[0], sample_regex)
        tool = normalize_tool(row.iloc[tool_col]) if tool_col is not None else default_tool
        cns = {}
        for locus, j in cn_cols.items():
            v = row.iloc[j]
            if not is_blank(v):
                try:
                    cns[canonical_locus(locus)] = int(float(v))
                except (TypeError, ValueError):
                    pass
        for locus, j in call_cols.items():
            raw = row.iloc[j]
            if is_blank(raw) and canonical_locus(locus) not in cns:
                continue
            out.extend(
                _make_records(sample, tool or "?", locus, raw, cns, sheet)
            )
    return out


def _read_matrix_layout(df, header_idx, sheet, default_tool, sample_regex,
                        cn_lookup=None, offset=0) -> list[Record]:
    header = [str(v) if pd.notna(v) else "" for v in df.iloc[header_idx]]
    loci = {
        j + offset: find_locus_in_text(h)
        for j, h in enumerate(header)
        if find_locus_in_text(h) and j + offset < df.shape[1]
    }
    out: list[Record] = []
    for i in range(header_idx + 1, len(df)):
        row = df.iloc[i]
        if is_blank(row.iloc[0]):
            continue
        sample = normalize_sample(row.iloc[0], sample_regex)
        cns = (cn_lookup or {}).get(sample, {})
        for j, locus in loci.items():
            raw = row.iloc[j]
            if is_blank(raw):
                continue
            out.extend(_make_records(sample, default_tool, locus, raw, cns, sheet))
    return out


def _read_wide_layout(df, header_idx, sheet, default_tool, sample_regex) -> list[Record]:
    """Layout 'planilha limpa': copy_number_X seguido de N colunas 'X alelo k'."""
    header = [str(v) if pd.notna(v) else "" for v in df.iloc[header_idx]]
    tool_col = next(
        (j for j, h in enumerate(header)
         if h.strip().lower() in {"software", "tool", "ping or kir-mapper"}),
        None,
    )

    # agrupa colunas por locus, na ordem em que aparecem
    blocks: list[tuple[str, int, list[int]]] = []  # (locus, cn_col, allele_cols)
    cur_locus, cur_cn, cur_alleles = None, None, []
    for j, h in enumerate(header):
        locus = find_locus_in_text(h)
        if locus is None:
            continue
        is_cn = "copy" in h.lower()
        if is_cn:
            if cur_locus:
                blocks.append((cur_locus, cur_cn, cur_alleles))
            cur_locus, cur_cn, cur_alleles = locus, j, []
        else:
            if cur_locus is None or canonical_locus(locus) != canonical_locus(cur_locus):
                if cur_locus:
                    blocks.append((cur_locus, cur_cn, cur_alleles))
                cur_locus, cur_cn, cur_alleles = locus, None, []
            cur_alleles.append(j)
    if cur_locus:
        blocks.append((cur_locus, cur_cn, cur_alleles))

    out: list[Record] = []
    for i in range(header_idx + 1, len(df)):
        row = df.iloc[i]
        if is_blank(row.iloc[0]):
            continue
        sample = normalize_sample(row.iloc[0], sample_regex)
        tool = normalize_tool(row.iloc[tool_col]) if tool_col is not None else default_tool
        for locus, cn_col, allele_cols in blocks:
            cn = None
            if cn_col is not None and not is_blank(row.iloc[cn_col]):
                try:
                    cn = int(float(row.iloc[cn_col]))
                except (TypeError, ValueError):
                    cn = None
            cells = [row.iloc[j] for j in allele_cols if not is_blank(row.iloc[j])]
            if not cells and cn is None:
                continue
            # cada celula = um alelo (possivelmente ambiguo com '/')
            raw = " + ".join(f"({str(c).strip()})" for c in cells) if cells else ""
            raw = raw.replace("(", "").replace(")", "")
            cns = {canonical_locus(locus): cn} if cn is not None else {}
            out.extend(_make_records(sample, tool or "?", locus, raw, cns, sheet))
    return out


def _read_cn_layout(df, header_idx, sample_regex, offset=0) -> dict[str, dict[str, int]]:
    header = [str(v) if pd.notna(v) else "" for v in df.iloc[header_idx]]
    loci = {
        j + offset: find_locus_in_text(h)
        for j, h in enumerate(header)
        if find_locus_in_text(h) and j + offset < df.shape[1]
    }
    out: dict[str, dict[str, int]] = {}
    for i in range(header_idx + 1, len(df)):
        row = df.iloc[i]
        if is_blank(row.iloc[0]):
            continue
        sample = normalize_sample(row.iloc[0], sample_regex)
        d = out.setdefault(sample, {})
        for j, locus in loci.items():
            v = row.iloc[j]
            if is_blank(v):
                continue
            try:
                d[canonical_locus(locus)] = int(float(v))
            except (TypeError, ValueError):
                pass
    return out


def _make_records(sample, tool, locus, raw, cns, sheet) -> list[Record]:
    """Cria um Record por gene canonico, separando loci compostos."""
    raw_s = "" if is_blank(raw) else str(raw).strip()
    if locus in COMPOSITE_LOCI and locus != "KIR2DL5":
        members = [canonical_locus(m) for m in COMPOSITE_LOCI[locus]]
        gts, warns = parse_call_string(raw_s, locus, copy_number=None)
        split = split_composite_locus(gts, locus, cns)
        recs = []
        for m in members:
            cn = cns.get(m)
            g = split.get(m, set())
            if cn is not None and g:
                from .parsing import _enforce_copy_number
                g = _enforce_copy_number(g, m, cn, warns)
            recs.append(Record(sample, tool, m, cn, raw_s, g, list(warns), sheet))
        return recs

    canon = canonical_locus(locus)
    cn = cns.get(canon)
    gts, warns = parse_call_string(raw_s, canon, copy_number=cn)
    return [Record(sample, tool, canon, cn, raw_s, gts, warns, sheet)]


# ------------------------------------------------------------------ API


def load_workbook(
    path: str | Path,
    sheet_tools: dict[str, str] | None = None,
    sample_regex: str = DEFAULT_SAMPLE_REGEX,
    skip_sheets: tuple[str, ...] = (),
) -> pd.DataFrame:
    """
    Le todas as abas relevantes e devolve a tabela longa de chamadas.

    `sheet_tools` mapeia nome-da-aba -> ferramenta, para abas que nao tem
    coluna de software (ex.: {'Sheet5': 'kir-mapper'}).
    """
    return load_sources([path], tool_map=sheet_tools, sample_regex=sample_regex,
                        skip=skip_sheets)


#: nomes de arquivo que denunciam a ferramenta de origem
_PADROES_FERRAMENTA = [
    ("PING", ("finalallelecalls", "iterallelecalls", "ping", "manualcopynumber",
              "locuscopynumber", "predictedcopynumber")),
    ("kir-mapper", ("kir-mapper", "kirmapper", "genotype", "ncopy", "calls")),
]

_EXTENSOES = {".csv", ".tsv", ".txt", ".xlsx", ".xls"}


def _ferramenta_pelo_nome(nome: str) -> str | None:
    n = nome.lower()
    if "ping" in n and "mapper" not in n:
        return "PING"
    if "mapper" in n and "ping" not in n:
        return "kir-mapper"
    for ferramenta, chaves in _PADROES_FERRAMENTA:
        if any(k in n for k in chaves):
            return ferramenta
    return None


def _quadros_do_arquivo(path: Path):
    """Devolve [(nome_do_bloco, DataFrame)] para xlsx (uma aba cada) ou csv/tsv."""
    if path.suffix.lower() in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(path)
        return [(sheet, xl.parse(sheet, header=None)) for sheet in xl.sheet_names]
    for sep in ("\t", ",", ";", r"\s+"):
        try:
            df = pd.read_csv(path, sep=sep, header=None, engine="python",
                             dtype=str, keep_default_na=False, na_values=[""])
        except Exception:
            continue
        if df.shape[1] > 1:
            return [(path.name, df)]
    return []


def load_sources(
    paths,
    tool_map: dict[str, str] | None = None,
    sample_regex: str = DEFAULT_SAMPLE_REGEX,
    skip: tuple[str, ...] = (),
) -> pd.DataFrame:
    """
    Le qualquer combinacao de arquivos e diretorios com saidas de PING e
    kir-mapper (xlsx, csv, tsv) e devolve a tabela longa de chamadas.

    A deteccao de layout e a mesma usada na planilha original, entao formatos
    novos nao exigem parser novo: basta que a tabela tenha uma coluna de amostra
    e colunas identificaveis por nome de locus KIR.

    `tool_map` associa nome de arquivo OU de aba a uma ferramenta. Quando
    ausente, a ferramenta e inferida do nome do arquivo/aba.
    """
    tool_map = tool_map or {}
    arquivos: list[Path] = []
    for entrada in ([paths] if isinstance(paths, (str, Path)) else paths):
        entrada = Path(entrada)
        if entrada.is_dir():
            arquivos += sorted(
                f for f in entrada.rglob("*") if f.suffix.lower() in _EXTENSOES
            )
        elif entrada.exists():
            arquivos.append(entrada)

    quadros: list[tuple] = []
    for arq in arquivos:
        for nome, df in _quadros_do_arquivo(arq):
            if nome in skip or df.empty:
                continue
            rotulo = nome if nome != arq.name else arq.name
            ferramenta = (tool_map.get(nome) or tool_map.get(arq.name)
                          or _ferramenta_pelo_nome(nome)
                          or _ferramenta_pelo_nome(arq.name)
                          or _ferramenta_pelo_nome(arq.parent.name))
            quadros.append((rotulo, df, ferramenta))

    return _processar_quadros(quadros, sample_regex)


def _processar_quadros(quadros, sample_regex: str) -> pd.DataFrame:
    """Duas passagens: primeiro as tabelas de copy number, depois as chamadas."""
    records: list[Record] = []
    cn_tables: dict[str, dict[str, dict[str, int]]] = {}
    pendentes: list[tuple] = []

    for nome, df, tool in quadros:
        hi = _header_row(df)
        if hi is None:
            continue
        layout = detect_layout(df, hi)
        off = header_offset(df, hi)
        if layout == "cn":
            alvo = tool or nome
            tabela = _read_cn_layout(df, hi, sample_regex, off)
            cn_tables.setdefault(alvo, {})
            for amostra, d in tabela.items():
                cn_tables[alvo].setdefault(amostra, {}).update(d)
        else:
            pendentes.append((nome, df, hi, layout, tool, off))

    for nome, df, hi, layout, tool, off in pendentes:
        cn_lookup = cn_tables.get(tool, {})
        if layout == "calls":
            recs = _read_calls_layout(df, hi, nome, tool, sample_regex)
        elif layout == "wide":
            recs = _read_wide_layout(df, hi, nome, tool, sample_regex)
        else:
            recs = _read_matrix_layout(df, hi, nome, tool, sample_regex, cn_lookup, off)
        records.extend(recs)

    return records_to_frame(records)


def _guess_tool_from_sheet(sheet: str) -> str | None:
    s = sheet.lower()
    if "ping" in s and "mapper" not in s:
        return "PING"
    if "mapper" in s and "ping" not in s:
        return "kir-mapper"
    return None


def records_to_frame(records: list[Record]) -> pd.DataFrame:
    rows = []
    for r in records:
        rows.append({
            "sample": r.sample,
            "tool": r.tool,
            "locus": r.locus,
            "copy_number": r.copy_number,
            "raw_call": r.raw_call,
            "n_candidates": len(r.genotypes),
            "genotypes": r.genotypes,
            "warnings": "; ".join(r.warnings),
            "source_sheet": r.source_sheet,
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["sample", "locus", "tool"]).reset_index(drop=True)
