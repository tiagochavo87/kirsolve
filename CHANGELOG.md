# Changelog

Registro das decisões de método e dos erros corrigidos durante o desenvolvimento.
Cada erro aqui alterou resultados, então o registro é parte da reprodutibilidade:
quem repetir a análise com uma versão anterior obtém números diferentes.

## [0.1.0] — 2026-08

### Adicionado

- Parsing de nomenclatura KIR (IPD-KIR), com loci compostos (`2DL23`, `2DS35`,
  `3DL1S1`, `2DL5AB`) e tokens `null` / `unresolved` / `failed`.
- Leitura de planilhas Excel multi-aba e de saídas nativas em CSV/TSV, com
  detecção automática de layout (`calls`, `matrix`, `cn`, `wide`).
- Consolidação intra-ferramenta e cruzamento PING × kir-mapper.
- Restrição por copy number.
- Resolução segura por nível de campo (resultado determinístico).
- Estimação de frequências alélicas por EM, com prior de Dirichlet.
- Análise de ubiquidade alélica e recomendação de nível de relato por locus.
- Validação: mascaramento leave-one-out, Hardy-Weinberg, comparação com
  frequências externas e comparação com verdade conhecida.
- Subcomando `benchmark` e conversor das anotações Immuannot (HPRC).

### Corrigido

**União trocada por interseção graduada no cruzamento entre ferramentas.**
A primeira versão, ao encontrar discordância entre PING e kir-mapper, unia os
candidatos das duas. Isso *aumentava* a ambiguidade em vez de reduzi-la: num
conjunto de 54 chamadas tipáveis, a soma de candidatos após o cruzamento era
393, contra 323 do PING sozinho e 154 do kir-mapper sozinho. O cruzamento
estava piorando o resultado.

A correção tenta a interseção em resolução menor (campo 2, depois campo 1)
antes de desistir. Discordância no terceiro campo costuma ser artefato de
resolução, não conflito real. Efeito: candidatos totais de 393 para 317,
ambiguidades de 21 para 7.

**Conflito irreconciliável deixou de ser rotulado como resolvido.**
A correção acima, sozinha, introduziu um erro pior: quando as ferramentas não
se reconciliavam nem em campo 1, o pipeline ficava com o conjunto menor e o
rotulava `resolvido_unico`. Os resolvidos saltaram de 26 para 36 — dez a mais,
todos falsos. Escolher o conjunto menor é parcimônia, não evidência.

Esses casos agora recebem o rótulo próprio `conflito_nao_resolvido` e exigem
inspeção do alinhamento.

**Deslocamento de coluna no cabeçalho.**
Abas sem a célula vazia no canto superior esquerdo têm os rótulos deslocados
uma coluna em relação aos dados, o que atribuía os alelos aos genes errados.
Detectado e corrigido automaticamente (`header_offset`).

**Ordem de alternância na regex de locus.**
`[23]D[LSP](?:\d[AB]?|\d\d|1S1)` casava `2DL2` dentro de `2DL23`, porque a
alternância é ordenada. Reordenado para `(?:1S1|\d[AB]{1,2}|\d\d|\d)`.

**Espaços ao redor de `+` e `/`.**
`"2DL3*001 / 2DL3*002 + 2DL3*008"` era quebrado pelo separador de espaço antes
de chegar ao parser de genótipo. Espaços adjacentes a `+` e `/` passaram a ser
normalizados antes da divisão.

**Critério de acerto na validação por mascaramento.**
O experimento comparava genótipos por igualdade literal, contando
`KIR2DL3*001` versus `KIR2DL3*0010101` como erro — são o mesmo alelo em
resoluções diferentes. Isso produziu acurácia de 8%, abaixo do acaso. Com
comparação por compatibilidade de prefixo, o valor real é 34,8% (campo 1) e
56,2% (campo 2), contra 25,3% e 38,2% de acerto ao acaso.

Verificado depois que o mesmo defeito não existia no pipeline principal: zero
pares redundantes em 80 conjuntos de candidatos, porque a consolidação passa
por `merge_genotypes` e a validação não passava.

**Classificação binária da verdade do HPRC.**
O conversor do Immuannot tratava qualquer `template_distance ≥ 1` como alelo
novo, o que classificava 86,5% dos loci como inalcançáveis. A distribuição real
é bimodal: um grupo de 0 a 3 bases (versão de banco ou variação não
codificante) e outro de 8 a 154 bases (alelo genuinamente novo). Limiar fixado
em 5, com três faixas em vez de duas.

### Limitações conhecidas

- O EM exige coorte grande. Abaixo de ~20 amostras informativas por locus as
  chamadas probabilísticas são suprimidas, e apenas o resultado determinístico
  vale.
- O prior de frequências precisa ser da população estudada. Frequências KIR
  variam muito entre grupos biogeográficos.
- Em três genomas de referência do HPRC, 54% dos loci KIR carregam alelo ausente
  do IPD-KIR. Nenhuma ferramenta baseada em catálogo pode acertá-los, e isso
  limita a base utilizável de qualquer benchmark.
- `2DL5A` e `2DL5B` são colapsados em `KIR2DL5`, porque nenhuma das duas
  ferramentas os separa de forma confiável em leitura curta.
- O modelo assume união aleatória de gametas. Coorte estruturada ou com
  aparentados enviesa as frequências.
