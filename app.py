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

# Em execucao local o pacote esta em src/; no navegador (stlite) ele vem
# instalado pelo wheel. O insert so acontece se a pasta existir, para nao
# mascarar o pacote instalado.
_src = Path(__file__).parent / "src"
if _src.is_dir():
    sys.path.insert(0, str(_src))

from kirsolve.loaders import load_sources
from kirsolve.report import write_outputs
from kirsolve.resolve import (
    consolidate,
    cross_tool,
    resolution_recommendation,
    resolve_cohort,
    ubiquity_report,
)

LOGO = Path(__file__).parent / "assets" / "logo.jpeg"

st.set_page_config(
    page_title="kirsolve",
    page_icon=str(LOGO) if LOGO.exists() else "🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

#: verde-petroleo do logo. Uma cor de acento so - o resto e neutro, para que
#: o destaque signifique alguma coisa quando aparecer.
ACENTO = "#3f7f6f"

st.markdown(f"""
<style>
  /* respiro no topo: o padrao do Streamlit cola o conteudo na barra */
  .block-container {{ padding-top: 2.5rem; max-width: 1200px; }}

  /* numeros em fonte tabular, para as colunas alinharem verticalmente */
  [data-testid="stMetricValue"] {{
      font-variant-numeric: tabular-nums;
      font-size: 1.9rem;
  }}
  [data-testid="stMetricLabel"] {{ color: #52525b; }}

  /* o cabecalho da secao precisa de hierarquia visivel sem virar enfeite */
  h3 {{ font-weight: 600; letter-spacing: -0.01em; }}

  /* area de upload: alvo grande e obvio, com a cor do acento na borda */
  [data-testid="stFileUploaderDropzone"] {{
      border: 2px dashed {ACENTO}55;
      background: {ACENTO}08;
      padding: 2rem 1rem;
  }}

  /* tabelas: numeros alinhados */
  [data-testid="stDataFrame"] {{ font-variant-numeric: tabular-nums; }}
</style>
""", unsafe_allow_html=True)

#: rotulo, o que fazer, e se exige acao do usuario.
#: A distincao util e binaria - da para usar, ou precisa de conferencia.
#: Uma escala de cinco cores obrigaria o leitor a decorar uma legenda antes
#: de ler o resultado, e o unico julgamento que ele precisa fazer e esse.
EXPLICACAO_DECISAO = {
    "resolvido_unico":        ("Resposta final",  "Alelo identificado com certeza total.", False),
    "resolvido_campo2":       ("Resposta parcial", "Certo até dois campos do nome. Não force mais detalhe.", False),
    "resolvido_campo1":       ("Resposta parcial", "Certo só na família do alelo.", False),
    "provavel_EM_alta":       ("Provável",        "Estimativa estatística, não observação. Confira antes de publicar.", True),
    "provavel_EM_moderada":   ("Hipótese",        "Estimativa estatística fraca. Trate como hipótese.", True),
    "absent":                 ("Gene ausente",    "A pessoa não tem esse gene. Isso é normal em KIR.", False),
    "no_call":                ("Sem dados",       "Nenhuma ferramenta conseguiu tipar.", False),
    "conflito_nao_resolvido": ("Conflito",        "As ferramentas discordam. Precisa conferência manual.", True),
    "ambiguo":                ("Sem resposta",    "Não foi possível decidir entre os candidatos.", True),
}


# ------------------------------------------------------------------ cabecalho

if LOGO.exists():
    st.image(str(LOGO), width=300)
st.markdown("#### Resolução de alelos KIR ambíguos")
st.caption(
    "Cruza as saídas do PING e do kir-mapper para reduzir a ambiguidade, "
    "e separa o que é certo do que é apenas provável."
)
st.divider()

# ------------------------------------------------------------------ barra lateral

with st.sidebar:
    if LOGO.exists():
        st.image(str(LOGO), use_container_width=True)
    st.caption("Sua planilha não sai deste computador. "
               "A análise roda inteiramente aqui.")
    st.divider()
    st.subheader("Opções")
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
    st.markdown("**Arraste a sua planilha para a área acima.**")
    st.caption("É o único passo necessário para começar.")
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

# Um diretorio por sessao, nao por rerun: o Streamlit reexecuta o script inteiro
# a cada interacao do usuario, e um mkdtemp() direto aqui acumularia uma pasta
# orfa por interacao. Guardado em session_state, o TemporaryDirectory e limpo
# sozinho quando a sessao termina (seu finalizador roda no garbage collection).
if "workdir" not in st.session_state:
    st.session_state.workdir = tempfile.TemporaryDirectory(prefix="kirsolve_")
pasta = Path(st.session_state.workdir.name)
for a in arquivos:
    (pasta / a.name).write_bytes(a.getbuffer())

# O erro aparece logo abaixo da area de envio, que e onde a acao aconteceu.
with st.spinner("Lendo os arquivos..."):
    try:
        bruto = load_sources([pasta], sample_regex=sample_regex)
    except Exception as exc:
        st.error(f"Não consegui ler os arquivos: {exc}")
        st.caption("Envie outro arquivo, ou confira se a planilha não está "
                   "aberta no Excel — isso pode bloquear a leitura.")
        st.stop()

if bruto.empty:
    st.error("Nenhuma chamada de alelo foi reconhecida nesses arquivos.")
    st.caption(
        "A planilha precisa de uma coluna com o nome da amostra e colunas "
        "nomeadas por gene KIR (KIR2DL1, KIR3DL3...). Se os nomes estiverem "
        "corretos, ajuste **Como encurtar o nome da amostra** na barra lateral."
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

pendentes = int(chamadas["decision"].isin(["ambiguo", "conflito_nao_resolvido"]).sum())
ausentes = int((chamadas["decision"] == "absent").sum())

m1, m2, m3 = st.columns([2, 2, 3])
m1.metric("Com resposta", resolvidos,
          delta=f"{resolvidos/n_tip:.0%} dos tipáveis" if n_tip else None,
          help="Alelo identificado pelo cruzamento entre as ferramentas.")
m2.metric("Precisam de atenção", pendentes,
          delta=None if not pendentes else "ver aba ao lado",
          delta_color="off",
          help="Ambíguos ou com conflito entre as ferramentas.")
m3.caption(
    f"De {len(chamadas)} análises, {ausentes} são genes que a pessoa não tem "
    f"— ausência de gene KIR é normal e não conta como falha. "
    f"Restam {n_tip} que precisavam de resposta."
)

if n_tip:
    st.progress(resolvidos / n_tip)

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
        {"decisão": k, "significa": rot,
         "exige conferência": "sim" if acao else "não", "o que fazer": desc}
        for k, (rot, desc, acao) in EXPLICACAO_DECISAO.items()
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
        column_config={
            "cópias": st.column_config.NumberColumn(format="%d", width="small"),
            "candidatos": st.column_config.NumberColumn(format="%d", width="small"),
            "probabilidade": st.column_config.NumberColumn(format="%.3f"),
        },
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
