#!/usr/bin/env bash
# =============================================================================
# Seleciona amostras que tenham SIMULTANEAMENTE:
#   - anotacao de verdade no HPRC (Zenodo 8372992)
#   - leitura curta aberta no 1000 Genomes
#   - a mesma populacao (exigencia do ncopy do kir-mapper)
#
# Uso: ./selecionar_amostras.sh [POPULACAO] [QUANTAS]
#      ./selecionar_amostras.sh YRI 50
#      ./selecionar_amostras.sh ""  50    # qualquer populacao, so para testar
#
# Por que uma populacao so: o autor do kir-mapper e explicito de que os
# limiares de copy number variam muito entre grupos biogeograficos, e
# recomenda rodar `map` e `ncopy` separados por ancestralidade. Misturar
# invalida a estimativa de numero de copias, que por sua vez determina o
# tamanho do genotipo. Ou seja: erra tudo o que vem depois.
# =============================================================================
set -euo pipefail

POP=${1:-}
QUANTAS=${2:-50}
BASE=${BASE:-/dados}

cd "$BASE"

# ---------------------------------------------------------------- verdade HPRC
if [[ ! -d verdade/hprc ]]; then
  echo "[1/4] baixando anotacoes do HPRC..."
  mkdir -p verdade && cd verdade
  curl -L -o hprc.tar "https://zenodo.org/records/8372992/files/hprc.tar?download=1"
  tar -xf hprc.tar
  cd "$BASE"
fi
ls verdade/hprc/*.gtf.gz | sed 's|.*/||; s|\.[12]\.gtf\.gz$||' | sort -u > /tmp/com_verdade.txt
echo "     $(wc -l < /tmp/com_verdade.txt) amostras com verdade"

# ---------------------------------------------------------------- indices 1000G
echo "[2/4] baixando indices do 1000 Genomes..."
for f in 1000G_2504_high_coverage 1000G_698_related_high_coverage; do
  [[ -f ${f}.sequence.index ]] || wget -q -c \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/${f}.sequence.index"
done

# painel de populacoes
[[ -f populacoes.txt ]] || wget -q -O populacoes.txt \
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped" \
  || echo "aviso: painel de populacoes indisponivel; filtro por populacao desativado" >&2

# ---------------------------------------------------------------- intersecao
echo "[3/4] cruzando..."
> /tmp/candidatas.txt
while read -r ID; do
  URL=$(grep -ho "ftp://[^[:space:]]*${ID}\.final\.cram" *.sequence.index 2>/dev/null | head -1 || true)
  [[ -z $URL ]] && continue
  if [[ -n $POP && -f populacoes.txt ]]; then
    LINHA=$(grep -P "\t${ID}\t" populacoes.txt 2>/dev/null | head -1 || true)
    [[ -z $LINHA ]] && continue
    # o .ped nao tem coluna de populacao em todas as versoes; se nao achar, mantem
    if echo "$LINHA" | grep -qv "$POP"; then
      echo "$LINHA" | grep -q "$POP" || continue
    fi
  fi
  printf '%s\t%s\n' "$ID" "$URL" >> /tmp/candidatas.txt
done < /tmp/com_verdade.txt

TOTAL=$(wc -l < /tmp/candidatas.txt)
echo "     $TOTAL amostras com verdade E leitura curta"

if (( TOTAL < QUANTAS )); then
  echo
  echo "AVISO: so ha $TOTAL amostras, menos que as $QUANTAS pedidas."
  echo "       Abaixo de 50, o ncopy do kir-mapper nao foi testado e o PING"
  echo "       quebra na analise de silhueta. Considere ampliar a populacao."
  echo
fi

head -n "$QUANTAS" /tmp/candidatas.txt > "$BASE/urls.txt"
cut -f1 "$BASE/urls.txt" > "$BASE/amostras.txt"

echo "[4/4] gravado:"
echo "     $BASE/amostras.txt  ($(wc -l < "$BASE/amostras.txt") IDs)"
echo "     $BASE/urls.txt      (ID + URL do CRAM)"
echo
head -3 "$BASE/urls.txt"
echo "..."
