# Validação do kirsolve

Registro dos resultados obtidos e do procedimento para reproduzi-los.
Data: agosto de 2026. Versão: 0.1.0.

## Resumo

Em 40 amostras do 1000 Genomes com montagem de leitura longa no HPRC como
padrão-ouro:

- **89,9% de acurácia** nos 90 loci em que o alelo verdadeiro consta do IPD-KIR
- **67,7%** sobre toda a verdade, incluindo alelos ausentes do catálogo
- **70,3% de acerto no desempate estatístico**, contra 32,5% ao acaso (229 casos)

O fator que mais limita o resultado não é o método: 72% dos loci do HPRC
carregam alelo que não existe no IPD-KIR, e nenhuma ferramenta baseada em
catálogo pode acertá-los.

---

## 1. O que foi validado

Duas coisas distintas, com métodos diferentes:

| O que | Como | Depende de |
|---|---|---|
| O desempate estatístico funciona? | Mascaramento (leave-one-out) | Só dos próprios dados |
| As chamadas estão corretas? | Comparação com montagem de leitura longa | Verdade externa (HPRC) |

A primeira mede a contribuição específica do kirsolve. A segunda mede o
resultado final, que depende também da qualidade do PING e do kir-mapper.

---

## 2. Resultados

### 2.1 Mascaramento — a tese central do método

Cada chamada já resolvida tem a resolução apagada artificialmente. As
frequências alélicas são estimadas com **todas as outras** amostras, e se
verifica se o modelo recupera a resposta certa.

| | 5 amostras | 36 amostras |
|---|---:|---:|
| Casos testados | 23 | **229** |
| Acurácia | 34,8% | **70,3%** |
| Acerto ao acaso | 25,3% | 32,5% |
| **Ganho sobre o acaso** | +9,5 pp | **+37,8 pp** |
| Verdade entre as 2 primeiras | 69,6% | 86,5% |

A coluna da esquerda vem das 5 amostras iniciais; a da direita, de 36 amostras
do 1000 Genomes.

**Interpretação.** Com 4 amostras de treino o modelo mal supera o chute. Com 34,
o ganho quadruplica. Isso confirma o mecanismo previsto: o EM desambigua porque
as amostras resolvidas em outras posições da coorte fornecem a informação para
desempatar as ambíguas. Sem coorte, não há informação, e o pipeline
corretamente suprime as chamadas probabilísticas.

O acerto é julgado por compatibilidade de prefixo, não igualdade literal:
`*001` e `*0010101` são o mesmo alelo em resoluções diferentes.

### 2.2 Comparação com a verdade do HPRC

Amostras com montagem fasada de leitura longa anotada pelo Immuannot
(Zenodo 8372992) **e** leitura curta aberta no 1000 Genomes.

Dois recortes, com 40 amostras extraídas:

| | Toda a verdade | Só onde o acerto é possível |
|---|---:|---:|
| Pares com verdade | 107 | 90 |
| Cobertura | 86,9% | 76,7% |
| **Acurácia** | 67,7% | **89,9%** |
| Verdade entre candidatos | 62,6% | 82,2% |
| Verdade descartada | 30 | 7 |

A coluna da direita usa `--apenas-exatos`, que mantém apenas os loci cujo alelo
verdadeiro consta do IPD-KIR. A diferença entre as duas colunas — 22 pontos
percentuais — é o custo imposto pelo catálogo, não pelo método.

Por locus, sobre a verdade completa:

| Locus | n | Acurácia |
|---|---:|---:|
| KIR2DS3 | 1 | 1,00 |
| KIR2DS4 | 10 | 0,89 |
| KIR2DL3 | 10 | 0,88 |
| KIR2DP1 | 10 | 0,78 |
| KIR3DL1 | 10 | 0,70 |
| KIR3DP1 | 10 | 0,67 |
| KIR2DS2, KIR2DS5, KIR3DS1 | 3 | 0,50 |
| KIR3DL3 | 10 | 0,50 |
| KIR3DL2 | 10 | 0,43 |
| KIR2DL5 | 3 | 0,33 |
| KIR2DS1 | 3 | 0,00 |

**Ressalva importante:** estes números saíram com o PING **incompleto**. Ele
parou antes de gerar o `finalAlleleCalls.csv`, e o pipeline usou o
`iterAlleleCalls.csv`, menos refinado. Com o PING completo, a tendência é
melhorar.

### 2.3 O teto imposto pelo catálogo

Nas 10 amostras analisadas, **72% dos loci têm alelo ausente do IPD-KIR**
(`template_distance > 5`, chegando a 154 bases de diferença). Nenhuma
ferramenta baseada em catálogo pode acertá-los.

Isso limita qualquer benchmark: conte com perder cerca de metade a dois terços
dos loci. Para 300 pares utilizáveis, planeje 35 a 40 amostras.

---

## 3. Os três achados que mais melhoraram o resultado

### 3.1 A região de extração estava errada

O intervalo usado inicialmente, `chr19:54.600.000–55.400.000`, foi escolhido por
estimativa. O `select_dna.bed` do próprio kir-mapper define
`chr19:54.114.000–54.905.000`.

Faltavam **486 mil bases** na ponta centromérica — exatamente onde ficam
KIR3DL3, KIR2DS2, KIR2DL2/3 e KIR2DL5B. O efeito na razão de cobertura do
KIR2DL5AB:

| Amostra | Antes | Depois |
|---|---:|---:|
| HG00673 | 0,041 | **0,429** |
| HG00733 | 0,069 | **0,483** |
| HG01175 | 0,095 | **0,553** |
| HG01243 | 0,142 | **0,943** |
| HG01258 | 0,135 | **1,143** |

O corte para 1 cópia é 0,25. Antes, **nenhuma** das 37 amostras chegava lá, e o
kir-mapper concluía que ninguém tinha o gene. A verdade do HPRC diz o contrário.

O BED também inclui duas janelas em chr6, na região HLA, que a versão inicial
ignorava.

**Lição:** as ferramentas trazem as coordenadas corretas nos próprios arquivos
de configuração. Consulte-as antes de definir a região.

### 3.2 O limiar de copy number do KIR3DP1

O `ncopy` não conseguiu calibrar os limiares e usou os **valores padrão** do
programa. Para o KIR3DP1, o corte para 3 cópias era 1,2.

| Amostra | Razão de cobertura | CN atribuído | CN verdadeiro |
|---|---:|---:|---:|
| HG00621 | 1,06 | 2 | 2 ✓ |
| HG00735 | 1,45 | 3 | 2 ✗ |
| HG01109 | 1,38 | 3 | 2 ✗ |

As duas erradas caíam logo acima do corte. Subindo o limiar de 1,2 para 1,6, as
duas voltaram a CN=2, e a acurácia do KIR3DP1 foi **de 0,00 para 0,40**.

Como o copy number determina o tamanho do genótipo, ele contamina tudo o que
vem depois. **Nenhuma chamada deve ser considerada definitiva antes da revisão
manual dos limiares.**

### 3.3 A fonte dos dados brutos

O FTP da ENA limita conexões e derruba transferências longas: 7 amostras em
8h34, com dezenas de retentativas. Os mesmos arquivos estão no bucket público
`1000genomes` da AWS, que não impõe esse limite:

```
https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/{ERR}/{AMOSTRA}.final.cram
```

Com a mesma máquina e o mesmo comando: **40 amostras em 9 minutos**. Prefira o
S3 ao FTP.

---

## 4. Limitações declaradas

**Leituras não mapeadas foram descartadas.** A extração pega a região KIR
(`chr19:54.6–55.4 Mb`) mais os contigs alternativos, sem as leituras não
mapeadas do resto do genoma. A primeira versão as incluía, mas 80% eram
descartadas depois por serem órfãs — o par ficara fora do recorte — e a
varredura dos 18 GB derrubava a conexão FTP. O benchmark, portanto, avalia o
insumo padrão das ferramentas, não robustez a leituras extraviadas.

**População misturada.** O autor do kir-mapper recomenda rodar `ncopy` por
grupo de ancestralidade, porque os limiares variam entre grupos
biogeográficos. Das 49 amostras com verdade no HPRC, apenas 40 têm leitura
curta aberta, o que não permitiu filtrar por população mantendo número.

**Versões de banco diferentes.** A anotação do HPRC usa IPD-KIR v2.12.0; o
kir-mapper compilado usa v2.15. Parte das discordâncias no terceiro campo pode
vir daí.

**n modesto.** 90 pares na comparação restrita, 229 no mascaramento. Suficiente
para detectar padrões por locus, insuficiente para intervalo de confiança
estreito.

**O PING não completou.** Parou antes de gerar o `finalAlleleCalls.csv`; foi
usado o `iterAlleleCalls.csv`, menos refinado. Os números tendem a melhorar com
a execução completa.

**Sete casos de `verdade_descartada` sem diagnóstico.** A investigação
identificou três causas — copy number errado, discordância real entre as
ferramentas e falha do PING em pseudogenes — mas os sete restantes não foram
analisados caso a caso.

---

## 5. Como reproduzir

### 5.1 Preparar

```bash
# amostras com verdade E leitura curta
./selecionar_amostras.sh "" 50

# pipeline completo (extração, kir-mapper, PING, benchmark)
nohup ./executar_tudo.sh > logs/pipeline.log 2>&1 &
```

Acompanhe com `cat ESTADO.txt`. O script tem retomada: se cair, rode de novo
que continua de onde parou.

### 5.2 Revisar os limiares de copy number

**Esta etapa é obrigatória e não pode ser automatizada.**

Para cada locus com acurácia baixa:

```bash
# 1. qual coluna é o locus
head -1 kirmapper/ncopy/ratio_values.txt | tr '\t' '\n' | cat -n | grep 3DP1

# 2. razões de cobertura das amostras
cut -f1,19 kirmapper/ncopy/ratio_values.txt

# 3. limiar atual (ordem: CN0:CN1:CN2:CN3)
grep KIR3DP1 kirmapper/ncopy/thresholds.txt

# 4. permissão (o Docker cria os arquivos como root)
docker run --rm -v /dados:/dados kir-mapper chown -R $(id -u):$(id -g) /dados/kirmapper

# 5. ajustar
python3 -c "
p='kirmapper/ncopy/thresholds.txt'
l=[('KIR3DP1:0.25:0.5:1.6:1.9' if x.startswith('KIR3DP1:') else x)
   for x in open(p).read().splitlines()]
open(p,'w').write('\n'.join(l)+'\n')"

# 6. refazer
docker run --rm -v /dados:/dados kir-mapper kir-mapper ncopy    -output /dados/kirmapper -threads 8
docker run --rm -v /dados:/dados kir-mapper kir-mapper genotype -output /dados/kirmapper -threads 8
```

Procure amostras cuja razão caia logo acima ou abaixo de um corte. São elas que
mudam de classe com um ajuste pequeno.

### 5.3 Avaliar

```bash
# verdade completa
python3 immuannot_para_verdade.py verdade/hprc -o verdade.tsv --amostras prontas.txt

# só onde o acerto é possível
python3 immuannot_para_verdade.py verdade/hprc -o verdade_exata.tsv \
    --amostras prontas.txt --apenas-exatos

python -m kirsolve benchmark ping/ kirmapper/genotype/cds \
    --truth verdade_exata.tsv -o benchmark/
```

Rode com as duas verdades e compare. A diferença entre elas é o custo imposto
pelo catálogo, não pelo método.

### 5.4 Ordem de leitura das métricas

| Métrica | O que significa |
|---|---|
| `verdade_descartada` | **Primeiro.** A resposta certa estava entre os candidatos e foi eliminada. Erro estrutural. |
| `sensibilidade_candidatos` | Teto da acurácia. Se está baixo, o problema é a montante do kirsolve. |
| `acuracia` | Acerto entre as que receberam chamada |
| `cobertura` | Proporção que recebeu chamada |

Acurácia sem olhar a sensibilidade engana: se a verdade não está entre os
candidatos, o kirsolve não tinha como acertar.

---

## 6. Próximos passos, em ordem de retorno

**1. Fazer o PING completar.** É a única etapa que não chegou ao fim, e afeta
metade do cruzamento. Investigar o erro em `logs/ping.log` e reexecutar.

**2. Diagnosticar os 7 casos de `verdade_descartada`.** Um a um, sem hipótese
prévia. Comparar `verdade`, `reportado`, `raw_PING` e `raw_kir_mapper` no
`benchmark_detalhado.tsv`. Só depois procurar padrão.

**3. Revisar os limiares dos loci ainda fracos.** KIR2DS1 (0,00) e KIR2DL5
(0,33). Mesmo procedimento da seção 3.2.

**4. Rodar `ncopy` por população.** Se conseguir 30 ou mais de um só grupo
biogeográfico, os limiares ficam mais estáveis e a comparação, mais limpa.

**5. Só então, ampliar a coorte.** O `cpc-p1.tar` do mesmo registro Zenodo tem
mais haplótipos anotados. Vale checar quantos têm leitura curta aberta.

### O que não fazer

Não ajuste limiares olhando a acurácia do benchmark. Isso é otimizar contra o
conjunto de teste, e o número resultante não vale nada. Ajuste olhando os
gráficos de razão de cobertura, que é o critério que o autor do kir-mapper
define — e só depois meça.

Não trate o resultado das etapas automáticas como final. O `ESTADO.txt` marca
explicitamente que os limiares precisam de revisão. Enquanto não forem
revisados, os copy numbers são provisórios e tudo abaixo deles também.

---

## 7. Origem dos dados

| Recurso | Fonte |
|---|---|
| Verdade (KIR anotado) | Zenodo 8372992, `hprc.tar` — Immuannot sobre montagens HPRC |
| Leitura curta | 1000 Genomes 30x, bucket S3 público `1000genomes` |
| Referência | GRCh38 full analysis set plus decoy HLA |
| kir-mapper | github.com/erickcastelli/kir-mapper (Docker) |
| PING | github.com/Hollenbach-lab/PING (Apptainer, Sylabs) |

Nenhum dado de paciente foi usado nesta validação.
