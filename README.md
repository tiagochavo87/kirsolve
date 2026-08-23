# kirsolve

Pipeline em Python para levantar e identificar os alelos KIR mais prováveis nos
casos em que o **PING** e o **kir-mapper** não chegaram a uma chamada única —
sobretudo os alelos ubíquos, aqueles que aparecem em quase todo conjunto
ambíguo porque não são discrimináveis com a cobertura disponível.

O pipeline **não substitui** as ferramentas: ele lê as saídas das duas, cruza as
evidências, aplica as restrições que cada uma sozinha não aplica, e só então usa
um modelo estatístico para desempatar o que sobrou.

---

## Lógica em cinco etapas

Cada etapa apenas **restringe** o conjunto de genótipos candidatos. Nada é
inventado, e o resultado determinístico é sempre separado do probabilístico.

| # | Etapa | O que faz |
|---|-------|-----------|
| 1 | **Consolidação** | Funde as várias abas/relatórios da mesma ferramenta. Interseção por compatibilidade de prefixo, mantendo sempre a maior resolução disponível (`*003` + `*00302` → `*00302`). |
| 2 | **Cruzamento PING × kir-mapper** | Interseção entre as duas ferramentas. Classifica em `concordant`, `one_tool`, `discordant`, `absent`, `no_call`. |
| 3 | **Copy number** | Cada genótipo é forçado ao CN observado. Genótipos com mais alelos que o CN são descartados ou fatiados; `*null` não conta como cópia. |
| 4 | **Resolução segura** | Maior nível de campo (1, 2 ou 3) em que **todos** os candidatos colapsam num único genótipo. É o resultado certo, sem modelo. |
| 5 | **EM + posterior** | Só para o que continuou ambíguo. Estima frequências alélicas na coorte por Expectation-Maximization e ranqueia os candidatos por probabilidade posterior. |

### O modelo estatístico (etapa 5)

Para um locus com copy number `c`, o genótipo é um multiconjunto de `c` alelos.
Sob união aleatória de gametas:

```
P(g) = c! / Π nₐ! · Π pₐ^nₐ
```

Cada amostra contribui com seu conjunto de candidatos `Cᵢ`. O EM alterna:

- **E**: `w(i,g) = P(g) / Σ_{g' ∈ Cᵢ} P(g')`
- **M**: `pₐ = (αₐ + Σᵢ Σ_g w(i,g)·nₐ(g)) / (Σα + total de cópias)`

`α` é um prior de Dirichlet. É ele que carrega a informação externa
(frequências publicadas, catálogo IPD-KIR). Isso importa porque **o EM só
desambigua quando há N**: amostras resolvidas em outras posições da coorte é que
dão a informação para desempatar as ambíguas. Com poucas dezenas de amostras o
prior domina, e o pipeline avisa e suprime as chamadas `provavel_EM`.

---

## Interface na web (GitHub Pages)

A interface está publicada em `https://SEU_USUARIO.github.io/kirsolve/` e roda
**inteiramente dentro do navegador**, via stlite (Streamlit compilado para
WebAssembly com Pyodide). Não há servidor: o GitHub Pages serve apenas arquivos
estáticos, e o Python executa na máquina de quem acessa.

**A planilha nunca é enviada a servidor nenhum.** Foi por isso que escolhemos
esta arquitetura em vez do Streamlit Community Cloud, onde os dados subiriam
para um servidor de terceiros.

⚠ Mesmo assim, para dado de paciente vale a cautela: a segurança do stlite em
contexto de dado sensível ainda não foi auditada de forma independente. Para
dado identificável, prefira rodar localmente.

Custo: cerca de 50 MB baixados na primeira visita (runtime Pyodide + pandas +
numpy). Depois disso funciona até offline.

### Como publicar

1. Settings → Pages → Source: **GitHub Actions**
2. `git push` na branch `main`

O workflow `.github/workflows/pages.yml` gera o wheel do kirsolve a cada push e
publica junto com a página. Isso garante que a interface no ar sempre usa o
código mais recente — sem esse passo, o site ficaria preso numa versão antiga.

---

## Interface gráfica local (para quem não usa linha de comando)

```bash
pip install -e . streamlit
streamlit run app.py
```

Abre no navegador. O usuário arrasta a planilha, confere o que foi lido, aperta
um botão e baixa o relatório em Excel — sem digitar comando nenhum.

A interface mostra as amostras, os genes e as ferramentas detectadas **antes**
de rodar, porque o erro mais comum é o nome da amostra não casar entre as duas
ferramentas, e isso aparece nessa conferência.

Toda a análise usa as mesmas funções da linha de comando (`load_sources`,
`consolidate`, `cross_tool`, `resolve_cohort`), então interface e CLI nunca
divergem.

Instruções passo a passo para o usuário final: `EXECUTAR_INTERFACE.txt`.

---

## Comece por aqui (sem precisar de dados)

O repositório inclui um gerador de dados sintéticos, então dá para rodar tudo
antes de ter qualquer arquivo real:

```bash
python scripts/gerar_exemplo.py exemplo --amostras 60
python -m kirsolve run exemplo/planilha_exemplo.xlsx -o resultados --sheet-tool "Sheet5=kir-mapper"
python -m kirsolve benchmark exemplo/ping exemplo/kirmapper \
    --truth exemplo/verdade.tsv -o benchmark --sample-regex '^(SIM\d+)'
```

Os dados são sorteados de um catálogo fixo com frequências arbitrárias. Não
representam população real e servem só para exercitar o código — mas cobrem os
casos que mais quebram parser: loci compostos, `unresolved`/`null`/`failed`,
anotação `$` de variante nova, cabeçalho deslocado e ambiguidade com `/`.

---

## Instalação

```bash
git clone <seu-repo> kirsolve && cd kirsolve
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e .
pytest -q          # 20 testes
```

### Docker

```bash
make docker
docker run --rm -v "$PWD":/work -w /work kirsolve:0.1.0 \
  run data/entrada.xlsx -o resultados --sheet-tool Sheet5=kir-mapper
```

### Servidor (SLURM)

```bash
sbatch scripts/slurm_kirsolve.sh dados.xlsx resultados/
```

O gargalo é memória, não CPU: a expansão combinatória dos candidatos é o que
pesa. Uma amostra com 19 alelos alternativos e CN=2 gera 190 genótipos; com CN=4
passaria de 8 mil. O parser tem teto configurável (`max_genotypes`).

---

## Uso

Antes de rodar, confira como o arquivo foi interpretado:

```bash
python -m kirsolve inspect dados.xlsx
```

Depois:

```bash
python -m kirsolve run dados.xlsx -o resultados/ \
    --sheet-tool "Sheet5=kir-mapper" \
    --priors frequencias_populacao.tsv \
    --prior-weight 5
```

| Opção | Para quê |
|---|---|
| `--sheet-tool ABA=FERRAMENTA` | Abas sem coluna de software. Repetível. |
| `--sample-regex` | Como extrair o ID curto. Padrão `^([A-Za-z]+\d+)`, que reduz `AMOSTRA001_23154FL-28Q2-01-140_S140_L004_R` → `AMOSTRA001`. |
| `--priors arquivo.tsv` | Frequências externas. Colunas: `locus`, `allele`, `frequency`. |
| `--prior-weight` | Peso do prior em "equivalente de amostras". |
| `--use-ipd` | Baixa o `Allelelist.txt` do IPD-KIR (ANHIG/IPDKIR) e sinaliza alelos fora do catálogo. |
| `--min-samples-em` | Abaixo disso o EM é marcado como não confiável. Padrão 20. |

### Entradas nativas das ferramentas

Além da planilha Excel, o pipeline lê diretórios com as saídas nativas:

```bash
python -m kirsolve inspect resultados_ping/
python -m kirsolve run resultados_ping/ resultados_kirmapper/ -o saida/
```

Não há parser específico por ferramenta. A mesma detecção automática de layout
serve para `finalAlleleCalls.csv` do PING, para a saída de `kir-mapper genotype`
e para formatos que ainda não vimos: basta a tabela ter uma coluna de amostra e
colunas identificáveis por nome de locus KIR.

⚠ Para benchmark, use o `finalAlleleCalls.csv` do PING, **não** o
`kirAllelesForAnalysis.csv`. O segundo já resolve ambiguidade por frequência
alélica — usá-lo mede o PING duas vezes e mascara o trabalho do kirsolve.

### Formatos de entrada reconhecidos

Detectados automaticamente, inclusive o deslocamento de coluna quando falta a
célula vazia do canto superior esquerdo:

- `Sample | software | <LOCUS>_Copy_number | <LOCUS>_Calls` (kir-mapper, `;`)
- `Sample | <LOCUS> | <LOCUS> | …` (PING, alternativas separadas por espaço)
- `Sample | <LOCUS> | …` só com inteiros → tabela de copy number
- `Sample | tool | copy_number_<X> | "<X> alelo 1" | "<X> alelo 2" | …`
  (planilha limpa, ambiguidade com `/`)

Loci compostos (`2DL23`, `2DS35`, `3DL1S1`, `2DL5AB`) são separados pelo próprio
nome do alelo. `2DL5A`/`2DL5B` são colapsados em `KIR2DL5`, porque nenhuma das
duas ferramentas os separa de forma confiável em leitura curta.

Anotações de variantes novas do PING (`KIR2DL3*0010101$E4_13.G^E4_15.A`) são
preservadas e vão para a aba de candidatos a alelo novo.

---

## Saídas

`kirsolve_relatorio.xlsx`, com as abas:

| Aba | Conteúdo |
|---|---|
| `00_resumo_por_locus` | Quanto o pipeline resolveu em cada locus |
| `01_chamadas_finais` | Uma linha por amostra × locus, com `decision`, `certain_call`, `best_call`, posterior e margem |
| `01b_matriz_chamadas` | Matriz amostra × locus, pronta para análise a jusante |
| `02_ainda_ambiguos` | O que sobrou sem resolver |
| `03_discordancias` | PING × kir-mapper sem nenhum genótipo compatível |
| `04_frequencias_EM` | Frequências estimadas, com flag de confiabilidade |
| `05_alelos_ubiquos` | Índice de ubiquidade por alelo |
| `06_resolucao_recomendada` | Em qual nível de campo vale reportar cada locus |
| `07_qc_copy_number` | CN divergente, gene framework com CN=0, CN atípico |
| `08_candidatos_novos` | Alelos com variantes extras anotadas pelo PING |
| `09_consolidado_por_ferramenta` | Rastreabilidade: chamada bruta, aba de origem, avisos |
| `10_diagnostico_EM` | Iterações, log-verossimilhança, convergência |

Os mesmos dados também saem em TSV.

### Como ler a coluna `decision`

| Valor | Significado |
|---|---|
| `resolvido_unico` | Um único genótipo em 3 campos. Confiança máxima. |
| `resolvido_campo2` / `campo1` | Certo até 2 ou 1 campo; ambíguo abaixo disso. **Reporte nesse nível**, não force resolução maior. |
| `provavel_EM_alta` | Posterior ≥ 0,90 e margem ≥ 0,50, com EM confiável. Chamada estatística, não observação. |
| `provavel_EM_moderada` | Posterior ≥ 0,70. Trate como hipótese. |
| `ambiguo` | Não resolvido. Precisa de mais evidência (cobertura, leitura longa, Sanger). |
| `absent` / `no_call` | CN=0 nas duas ferramentas / nenhuma tipou. |

### O índice de ubiquidade

```
ubiquity_index = prop_sets × (1 − n_certain / n_sets)
```

Próximo de 1: o alelo aparece em quase todo conjunto ambíguo mas quase nunca é
certo — é ruído de nomenclatura, um nome que a resolução disponível não
distingue dos vizinhos. Próximo de 0: alelo informativo.

Alelos com índice alto num locus são o argumento concreto para reportar aquele
locus em 2 campos em vez de 3. Use junto com a aba `06_resolucao_recomendada`.

---

## Validação

Quatro estratégias, em `kirsolve/validate.py`.

### Sem dados externos

**Mascaramento (validação cruzada leave-one-out).** Pega as chamadas já
resolvidas, apaga artificialmente a resolução, estima frequências com todas as
outras amostras e verifica se o modelo recupera a resposta certa. É a medida
direta da acurácia do desempate estatístico.

O acerto é julgado por **compatibilidade de prefixo**, não igualdade literal:
`*001` e `*0010101` são o mesmo alelo em resoluções diferentes. Medir por
igualdade produziria acurácia abaixo do acaso (ver `CHANGELOG.md`).

**Hardy-Weinberg.** Qui-quadrado por locus, usando amostras com CN=2 e genótipo
único. Exige ao menos 30 amostras. Desvio forte indica erro sistemático — alelo
nulo não detectado, ou leituras de parálogos contaminando a chamada.

**Comparação de frequências.** Confronta as estimativas com uma referência
externa. Divergência grande em alelos comuns é sinal de viés.

### Contra verdade conhecida

```bash
python -m kirsolve benchmark saidas_ping/ saidas_kirmapper/ \
    --truth verdade.tsv -o benchmark/
```

Saídas: `benchmark_resumo.json`, `benchmark_por_locus.tsv`,
`benchmark_detalhado.tsv`.

| Métrica | Leitura |
|---|---|
| `verdade_descartada` | **Olhe primeiro.** Se > 0, alguma etapa eliminou a resposta certa — erro estrutural, pior que errar a escolha. |
| `sensibilidade_candidatos` | Proporção em que a verdade estava entre os candidatos. É o teto da acurácia. |
| `acuracia` | Acerto entre as que receberam chamada |
| `cobertura` | Proporção que recebeu chamada |

A comparação respeita a resolução: uma chamada reportada em campo 2 é julgada
em campo 2. Cobrar três campos de quem só afirmou dois mediria a resolução do
sequenciamento, não a acurácia do método.

### Obtendo a verdade (HPRC)

As anotações do Immuannot para os genomas montados do HPRC estão públicas no
Zenodo, registro `8372992`, arquivo `hprc.tar` — um GTF por haplótipo. Não é
preciso rodar o Immuannot.

```bash
curl -L -o hprc.tar "https://zenodo.org/records/8372992/files/hprc.tar?download=1"
tar -xf hprc.tar
python scripts/immuannot_para_verdade.py hprc/ -o verdade.tsv --amostras amostras.txt
```

Como a montagem é fasada, cada amostra tem dois GTF (`.1` e `.2`) e o conversor
junta os dois num genótipo — a verdade já vem resolvida.

O conversor classifica cada locus em três faixas pelo `template_distance`, que
mede quantas bases a montagem difere do alelo de referência mais próximo:

| Faixa | Significado |
|---|---|
| 0 | idêntico ao IPD-KIR |
| 1–5 | versão de banco ou variação não codificante; não impede o acerto |
| > 5 | **alelo novo** — nenhuma ferramenta baseada em catálogo pode acertar |

O limiar foi calibrado na distribuição real, que é bimodal com vão entre 3 e 8.
Um corte binário em ≥ 1 classificaria 86% dos loci como inalcançáveis.

Em três genomas de referência (HG00733, NA19240, HG02818), **54% dos loci KIR
carregam alelo ausente do IPD-KIR** — um `KIR3DL2` chegou a 154 bases de
diferença. Isso limita a base utilizável de qualquer benchmark: conte com
perder cerca de metade dos loci. Para 300 pares utilizáveis, planeje 35–40
amostras.

---

## Limitações

- **O EM precisa de N.** Com cinco amostras ele não desambigua nada; só o prior
  externo trabalha. A partir de ~150–200 amostras da mesma população as
  frequências ficam úteis. O pipeline avisa e suprime chamadas EM abaixo do
  limiar.
- **O prior tem que ser da população certa.** Frequências KIR variam muito entre
  grupos biogeográficos. Um prior europeu numa coorte brasileira miscigenada
  enviesa a chamada. Sem prior adequado, prefira reportar em resolução menor.
- **Discordância não é erro do pipeline.** Quando PING e kir-mapper não têm
  genótipo compatível, o caso vai para `03_discordancias` com as chamadas brutas
  das duas ferramentas. Esses casos pedem inspeção do BAM — o kir-mapper permite
  isso e é o desempate mais confiável.
- **Alelos ubíquos podem ser irredutíveis.** Se dois alelos diferem só em região
  não coberta pelo sequenciamento, nenhum algoritmo resolve. A resposta correta
  é reportar em campo 1 ou 2, não escolher um.
- **HWE é uma aproximação.** O modelo assume união aleatória de gametas. Em
  coorte estruturada ou com aparentados, as frequências saem enviesadas.

## Estrutura do repositório

```
src/kirsolve/
  nomenclature.py  parsing de alelos, loci compostos, níveis de resolução
  parsing.py       strings de chamada das duas ferramentas
  loaders.py       leitura de xlsx/csv/tsv com detecção de layout
  resolve.py       consolidação, cruzamento, resolução segura
  em.py            Expectation-Maximization e posterior
  priors.py        frequências externas e catálogo IPD-KIR
  validate.py      mascaramento, HWE, comparação com verdade
  report.py        relatórios
  cli.py           run / inspect / benchmark
scripts/
  gerar_exemplo.py            dados sintéticos
  immuannot_para_verdade.py   GTF do HPRC -> tabela de verdade
  slurm_kirsolve.sh           submissão em cluster
tests/                        20 testes
```

`CHANGELOG.md` registra os erros de método corrigidos durante o
desenvolvimento, incluindo dois que alteraram resultados. Vale ler antes de
comparar números entre versões.

---

## Referências

- Castelli EC et al. *kir-mapper: A Toolkit for KIR Genotyping From Short-Read
  Second-Generation Sequencing Data.* HLA, 2025. doi:10.1111/tan.70092
- Norman PJ et al. *Defining KIR and HLA Class I Genotypes at Highest Resolution
  via High-Throughput Sequencing.* Am J Hum Genet, 2016.
- IPD-KIR Database — nomenclatura e formato GL string:
  https://www.ebi.ac.uk/ipd/kir/
- Catálogo de alelos: https://github.com/ANHIG/IPDKIR
