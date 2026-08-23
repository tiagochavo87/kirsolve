"""
Interface grafica do kirsolve, para quem nao usa linha de comando.

Como rodar:
    pip install streamlit
    streamlit run app.py

Abre no navegador. O usuario arrasta a planilha, confere o que foi lido,
aperta um botao e baixa o relatorio.

Toda a analise usa as mesmas funcoes do pacote (`load_sources`, `consolidate`,
`cross_tool`, `resolve_cohort`, `write_outputs`). Este arquivo e apenas a
camada visual - nao ha logica de analise duplicada aqui, para que a interface
e a linha de comando nunca divirjam.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "src"))

from kirsolve.loaders import load_sources
from kirsolve.report import write_outputs
from kirsolve.resolve import (
    consolidate,
    cross_tool,
    resolution_recommendation,
    resolve_cohort,
    ubiquity_report,
)

st.set_page_config(page_title="kirsolve", page_icon="🧬", layout="wide")

LOGO = Path(__file__).parent / "assets" / "logo.jpeg"

EXPLICACAO_DECISAO = {
    "resolvido_unico": ("Resposta final", "Alelo identificado com certeza total.", "🟢"),
    "resolvido_campo2": ("Resposta parcial", "Certo até dois campos do nome. Não force mais detalhe.", "🟢"),
    "resolvido_campo1": ("Resposta parcial", "Certo só na família do alelo.", "🟡"),
    "provavel_EM_alta": ("Provável", "Estimativa estatística com alta confiança. Não é observação.", "🟡"),
    "provavel_EM_moderada": ("Hipótese", "Estimativa estatística moderada. Trate como hipótese.", "🟠"),
    "absent": ("Gene ausente", "A pessoa não tem esse gene. Isso é normal em KIR.", "⚪"),
    "no_call": ("Sem dados", "Nenhuma ferramenta conseguiu tipar.", "⚪"),
    "conflito_nao_resolvido": ("Conflito", "As ferramentas discordam. Precisa conferência manual.", "🔴"),
    "ambiguo": ("Sem resposta", "Não foi possível decidir entre os candidatos.", "🔴"),
}


# ------------------------------------------------------------------ cabecalho

col_logo, col_txt = st.columns([1, 3])
with col_logo:
    if LOGO.exists():
        st.image(str(LOGO), width=260)
with col_txt:
    st.markdown("### Resolução de alelos KIR ambíguos")
    st.caption("Cruza as saídas do PING e do kir-mapper para reduzir a ambiguidade "
               "e separar o que é certo do que é apenas provável.")

st.divider()

# ------------------------------------------------------------------ barra lateral

with st.sidebar:
    st.header("Opções")
    st.caption("Os padrões funcionam na maioria dos casos.")

    sample_regex = st.text_input(
        "Como encurtar o nome da amostra",
        value=r"^([A-Za-z]+\d+)",
        help="Expressão regular. O padrão transforma "
             "'AMOSTRA001_23154FL-28Q2' em 'AMOSTRA001'. "
             "Serve para casar os nomes entre as duas ferramentas.",
    )

    min_em = st.number_input(
        "Mínimo de amostras para confiar na estatística",
        min_value=5, max_value=200, value=20, step=5,
        help="Abaixo disso, o desempate estatístico é suprimido e só valem "
             "as resoluções certas. Não recomendamos reduzir.",
    )

    st.divider()
    st.markdown("**Precisa de ajuda?**")
    st.caption("Se o resultado vier vazio, quase sempre é o nome da amostra "
               "que não casou entre as ferramentas. Confira a aba "
               "'O que foi lido' antes de rodar.")


# ------------------------------------------------------------------ upload

st.subheader("1. Envie os arquivos")

arquivos = st.file_uploader(
    "Planilha Excel, ou os arquivos CSV/TSV das ferramentas",
    type=["xlsx", "xls", "csv", "tsv", "txt"],
    accept_multiple_files=True,
    help="Pode enviar a planilha única com todas as abas, ou os arquivos "
         "separados do PING e do kir-mapper de uma vez.",
)

if not arquivos:
    st.info("Envie ao menos um arquivo para começar.")
    with st.expander("Que arquivos o programa aceita?"):
        st.markdown("""
- **Planilha Excel** com uma ou mais abas, uma coluna por gene
- **`finalAlleleCalls.csv`** do PING
- **`manualCopyNumberFrame.csv`** do PING
- **arquivos `.calls.txt`** do kir-mapper (um por gene)
- qualquer CSV ou TSV com uma coluna de amostra e colunas nomeadas por gene KIR

O programa reconhece o formato sozinho. Não é preciso renomear nada.
        """)
    st.stop()


# ------------------------------------------------------------------ leitura

pasta = Path(tempfile.mkdtemp(prefix="kirsolve_"))
for a in arquivos:
    (pasta / a.name).write_bytes(a.getbuffer())

with st.spinner("Lendo os arquivos..."):
    try:
        bruto = load_sources([pasta], sample_regex=sample_regex)
    except Exception as exc:
        st.error(f"Não consegui ler os arquivos: {exc}")
        st.stop()

if bruto.empty:
    st.error(
        "Nenhuma chamada de alelo foi reconhecida nesses arquivos.\n\n"
        "Causas comuns: a planilha não tem colunas com nome de gene KIR, "
        "ou a coluna de amostra está em outro lugar."
    )
    st.stop()

st.subheader("2. Confira o que foi lido")
st.caption("Antes de rodar, verifique se as amostras, os genes e as ferramentas "
           "estão como você espera.")

c1, c2, c3 = st.columns(3)
c1.metric("Amostras", bruto["sample"].nunique())
c2.metric("Genes", bruto["locus"].nunique())
c3.metric("Ferramentas", bruto["tool"].dropna().nunique())

ferramentas = sorted(bruto["tool"].dropna().unique())
if len(ferramentas) < 2:
    st.warning(
        f"Só encontrei uma ferramenta: **{ferramentas[0] if ferramentas else '?'}**. "
        "O cruzamento entre PING e kir-mapper é o que mais reduz ambiguidade — "
        "sem ele o resultado fica bem mais fraco. Confira se faltou enviar algum arquivo."
    )

with st.expander("Ver detalhes do que foi lido"):
    st.dataframe(
        bruto.groupby(["source_sheet", "tool"]).size().reset_index(name="registros"),
        use_container_width=True, hide_index=True,
    )
    st.write("**Amostras:**", ", ".join(sorted(bruto["sample"].unique())))
    st.write("**Genes:**", ", ".join(sorted(bruto["locus"].unique())))
    avisos = bruto[bruto["warnings"] != ""]
    if len(avisos):
        st.warning(f"{len(avisos)} linhas com avisos de leitura:")
        st.dataframe(avisos[["sample", "tool", "locus", "warnings"]].head(20),
                     use_container_width=True, hide_index=True)


# ------------------------------------------------------------------ analise

st.subheader("3. Analisar")

if not st.button("Rodar análise", type="primary", use_container_width=True):
    st.stop()

with st.spinner("Cruzando as ferramentas e resolvendo..."):
    consolidado = consolidate(bruto)
    cruzado = cross_tool(consolidado)
    chamadas, freqs, diags = resolve_cohort(cruzado, min_samples_for_em=int(min_em))
    ubiq = ubiquity_report(cruzado)
    resol = resolution_recommendation(cruzado)

    saida = pasta / "resultados"
    arquivos_gerados = write_outputs(
        saida, chamadas, freqs, ubiq, resol, cruzado, consolidado, diags
    )

st.success("Análise concluída.")

# ------------------------------------------------------------------ resultados

tipaveis = chamadas[~chamadas["decision"].isin(["absent", "no_call"])]
resolvidos = int(tipaveis["decision"].str.startswith("resolvido").sum())
n_tip = len(tipaveis)

m1, m2, m3, m4 = st.columns(4)
m1.metric("Análises", len(chamadas))
m2.metric("Genes ausentes", int((chamadas["decision"] == "absent").sum()),
          help="A pessoa não tem o gene. Normal em KIR.")
m3.metric("Com resposta", resolvidos,
          delta=f"{resolvidos/n_tip:.0%} dos tipáveis" if n_tip else None)
m4.metric("Precisam de atenção",
          int(chamadas["decision"].isin(["ambiguo", "conflito_nao_resolvido"]).sum()))

if not any(d.get("em_reliable") for d in diags.values()):
    st.info(
        f"**O desempate estatístico não foi usado.** Ele precisa de pelo menos "
        f"{int(min_em)} amostras para funcionar: são as amostras já resolvidas "
        "que dão a informação para desempatar as ambíguas. Todas as respostas "
        "aqui vêm do cruzamento entre as ferramentas, sem modelo estatístico."
    )

aba1, aba2, aba3, aba4 = st.tabs(
    ["Resultados", "Precisam de atenção", "Por gene", "Baixar"]
)

with aba1:
    st.caption("Uma linha por amostra e gene. A coluna **decisão** diz o quanto "
               "confiar; a coluna **alelo** traz a resposta.")

    legenda = pd.DataFrame([
        {"": ico, "decisão": k, "significa": rot, "o que fazer": desc}
        for k, (rot, desc, ico) in EXPLICACAO_DECISAO.items()
        if k in set(chamadas["decision"])
    ])
    with st.expander("O que cada decisão significa"):
        st.dataframe(legenda, use_container_width=True, hide_index=True)

    filtro = st.multiselect(
        "Mostrar apenas:", sorted(chamadas["decision"].unique()), default=[]
    )
    tabela = chamadas[chamadas["decision"].isin(filtro)] if filtro else chamadas

    st.dataframe(
        tabela[["sample", "locus", "copy_number", "decision",
                "certain_call", "posterior", "n_candidates"]]
        .rename(columns={
            "sample": "amostra", "locus": "gene", "copy_number": "cópias",
            "decision": "decisão", "certain_call": "alelo",
            "posterior": "probabilidade", "n_candidates": "candidatos",
        }),
        use_container_width=True, hide_index=True, height=420,
    )

with aba2:
    problemas = chamadas[chamadas["decision"].isin(["ambiguo", "conflito_nao_resolvido"])]
    if problemas.empty:
        st.success("Nenhum caso pendente.")
    else:
        st.caption("Nestes casos o programa não conseguiu decidir. "
                   "As colunas da direita mostram o que cada ferramenta disse.")
        st.dataframe(
            problemas[["sample", "locus", "decision", "n_candidates",
                       "raw_PING", "raw_kir_mapper"]]
            .rename(columns={
                "sample": "amostra", "locus": "gene", "decision": "motivo",
                "n_candidates": "candidatos",
                "raw_PING": "o que o PING disse",
                "raw_kir_mapper": "o que o kir-mapper disse",
            }),
            use_container_width=True, hide_index=True, height=420,
        )

with aba3:
    st.caption("Quanto o programa resolveu em cada gene, e em que nível de "
               "detalhe vale relatar.")
    from kirsolve.report import summarize
    resumo = summarize(chamadas)
    st.dataframe(resumo, use_container_width=True, hide_index=True)
    if not resol.empty:
        st.markdown("**Nível de detalhe recomendado por gene**")
        st.caption("Insistir em mais detalhe do que o recomendado introduz erro.")
        st.dataframe(resol, use_container_width=True, hide_index=True)

with aba4:
    st.caption("O arquivo Excel traz tudo, em abas separadas.")
    for f in arquivos_gerados:
        f = Path(f)
        st.download_button(
            f"Baixar {f.name}",
            data=f.read_bytes(),
            file_name=f.name,
            use_container_width=True,
        )

st.divider()
st.caption(
    "kirsolve · o resultado marcado como **resolvido** vem do cruzamento entre "
    "as ferramentas, não de estimativa. O marcado como **provável** vem de "
    "modelo estatístico e deve ser tratado como hipótese."
)
