#!/usr/bin/env bash
# =============================================================================
# Pipeline completo, para rodar a noite sem supervisao.
#
#   extracao -> kir-mapper (map, ncopy, genotype) -> PING -> kirsolve benchmark
#
# Uso:
#     nohup /dados/executar_tudo.sh > /dados/logs/pipeline.log 2>&1 &
#     tail -f /dados/logs/pipeline.log
#
# Estado, para acompanhar de fora:
#     cat /dados/ESTADO.txt
#
# -----------------------------------------------------------------------------
# DECISOES QUE IMPORTAM
#
# Retentativa: das tres amostras do teste piloto, DUAS falharam por queda do
# FTP da ENA no meio da transferencia. Sem retentativa automatica, uma queda
# as 3h da manha perde a noite inteira. Cada amostra tem 4 tentativas com
# espera crescente.
#
# Paralelismo limitado: puxar muitas amostras ao mesmo tempo do mesmo servidor
# da ENA aumenta a taxa de recusa - foi o que derrubou a HG02818 enquanto
# outro download corria. PARALELO=3 e um meio-termo entre velocidade e
# educacao com o servidor alheio.
#
# Retomada: tudo verifica se ja esta feito antes de refazer. Se cair, e so
# rodar de novo que continua de onde parou.
#
# Etapas manuais: `ncopy` (kir-mapper) e o copy number do PING pedem revisao
# humana dos limiares. O script roda a versao automatica e SEGUE, deixando
# marcado no ESTADO.txt que os limiares precisam de conferencia. Nao trave a
# noite esperando por isso, mas nao considere os copy numbers definitivos.
# =============================================================================
set -uo pipefail

BASE=${BASE:-/dados}
THREADS=${THREADS:-8}
PARALELO=${PARALELO:-3}
TENTATIVAS=${TENTATIVAS:-4}
REF="$BASE/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REGIAO="chr19:54600000-55400000"
ESTADO="$BASE/ESTADO.txt"

export PATH="$BASE/tools/bin:$PATH"
mkdir -p "$BASE"/{bam,fastq,kirmapper,ping,kirsolve,logs}

# ---------------------------------------------------------------- utilidades
marcar() { printf '%s | %s\n' "$(date '+%d/%m %H:%M')" "$*" | tee -a "$ESTADO"; }
titulo() { echo; echo "==================== $* ===================="; }

verificar_requisitos() {
  local faltando=0
  for cmd in samtools docker apptainer python3; do
    command -v $cmd >/dev/null || { echo "FALTA: $cmd"; faltando=1; }
  done
  [[ -f $REF ]] || { echo "FALTA: referencia em $REF"; faltando=1; }
  [[ -f $REF.fai ]] || { echo "FALTA: indice $REF.fai"; faltando=1; }
  [[ -f $BASE/urls.txt ]] || { echo "FALTA: $BASE/urls.txt (rode selecionar_amostras.sh)"; faltando=1; }
  docker images 2>/dev/null | grep -q kir-mapper || { echo "FALTA: imagem docker kir-mapper"; faltando=1; }
  [[ -f $BASE/tools/PING/ping.sif ]] || { echo "FALTA: ping.sif"; faltando=1; }
  return $faltando
}

# ---------------------------------------------------------------- 1. extracao
extrair_uma() {
  local ID=$1 URL=$2
  [[ -f $BASE/fastq/${ID}_R1.fq.gz ]] && return 0

  local ALT
  ALT=$(grep -E 'KI270921|KI270923|KI270890' "$REF.fai" | cut -f1 | tr '\n' ' ')

  local n=1
  while (( n <= TENTATIVAS )); do
    if samtools view -@ 2 -T "$REF" -b -o "$BASE/bam/${ID}.a.bam" "$URL" $REGIAO $ALT \
         2>> "$BASE/logs/extrair_${ID}.log" \
       && samtools view -@ 2 -T "$REF" -b -f 4 -o "$BASE/bam/${ID}.b.bam" "$URL" \
         2>> "$BASE/logs/extrair_${ID}.log"
    then
      samtools merge -@ 2 -f "$BASE/bam/${ID}.kir.bam" \
        "$BASE/bam/${ID}".{a,b}.bam 2>> "$BASE/logs/extrair_${ID}.log" || { n=$((n+1)); continue; }
      samtools sort -@ 2 -o "$BASE/bam/${ID}.sorted.bam" "$BASE/bam/${ID}.kir.bam" \
        2>> "$BASE/logs/extrair_${ID}.log"
      samtools index "$BASE/bam/${ID}.sorted.bam"
      samtools collate -@ 2 -u -O "$BASE/bam/${ID}.kir.bam" 2>/dev/null \
        | samtools fastq -@ 2 -1 "$BASE/fastq/${ID}_R1.fq.gz" \
                         -2 "$BASE/fastq/${ID}_R2.fq.gz" \
                         -0 /dev/null -s /dev/null -n 2>> "$BASE/logs/extrair_${ID}.log"
      rm -f "$BASE/bam/${ID}".{a,b,kir}.bam
      [[ -s $BASE/fastq/${ID}_R1.fq.gz ]] && { echo "  ok   $ID"; return 0; }
    fi
    rm -f "$BASE/bam/${ID}".{a,b,kir}.bam
    echo "  retry $ID (tentativa $n de $TENTATIVAS)"
    sleep $(( n * 60 ))
    n=$((n+1))
  done
  echo "  FALHOU $ID apos $TENTATIVAS tentativas"
  echo "$ID" >> "$BASE/logs/amostras_falhadas.txt"
  return 1
}
export -f extrair_uma
export BASE THREADS TENTATIVAS REF REGIAO

etapa_extracao() {
  titulo "1. EXTRACAO"
  local total pendentes
  total=$(wc -l < "$BASE/urls.txt")
  pendentes=$(while IFS=$'\t' read -r ID _; do
                [[ -f $BASE/fastq/${ID}_R1.fq.gz ]] || echo "$ID"
              done < "$BASE/urls.txt" | wc -l)
  marcar "extracao: $pendentes pendentes de $total"

  # xargs cuida do paralelismo; cada amostra tem retentativa propria
  while IFS=$'\t' read -r ID URL; do
    printf '%s\t%s\n' "$ID" "$URL"
  done < "$BASE/urls.txt" \
    | xargs -P "$PARALELO" -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r i u <<< "{}"; extrair_uma "$i" "$u"'

  local prontas
  prontas=$(ls "$BASE"/fastq/*_R1.fq.gz 2>/dev/null | wc -l)
  marcar "extracao concluida: $prontas de $total amostras"
  (( prontas > 0 ))
}

# ---------------------------------------------------------------- 2. kir-mapper
etapa_kirmapper() {
  titulo "2. KIR-MAPPER"
  local n=0
  for BAMF in "$BASE"/bam/*.sorted.bam; do
    [[ -e $BAMF ]] || continue
    local ID; ID=$(basename "$BAMF" .sorted.bam)
    [[ -d $BASE/kirmapper/map/$ID ]] && continue
    echo "  map $ID"
    docker run --rm -v "$BASE":"$BASE" kir-mapper kir-mapper map \
      -bam "$BAMF" -sample "$ID" -output "$BASE/kirmapper" -threads "$THREADS" \
      >> "$BASE/logs/map.log" 2>&1 && n=$((n+1))
  done
  marcar "kir-mapper map: $n novas amostras"

  echo "  ncopy"
  docker run --rm -v "$BASE":"$BASE" kir-mapper kir-mapper ncopy \
    -output "$BASE/kirmapper" -threads "$THREADS" >> "$BASE/logs/ncopy.log" 2>&1
  if [[ -f $BASE/kirmapper/ncopy/copy_numbers.table.txt ]]; then
    marcar "ncopy ok - LIMIARES PRECISAM DE REVISAO MANUAL (ver ncopy/plots/*.html)"
  else
    marcar "ncopy FALHOU - ver logs/ncopy.log"
    return 1
  fi

  echo "  genotype"
  docker run --rm -v "$BASE":"$BASE" kir-mapper kir-mapper genotype \
    -output "$BASE/kirmapper" -threads "$THREADS" >> "$BASE/logs/genotype.log" 2>&1
  local nc
  nc=$(find "$BASE/kirmapper/genotype" -name "*.calls.txt" 2>/dev/null | wc -l)
  marcar "kir-mapper genotype: $nc arquivos de chamada"
  (( nc > 0 ))
}

# ---------------------------------------------------------------- 3. PING
etapa_ping() {
  titulo "3. PING"
  if [[ -s $BASE/ping/finalAlleleCalls.csv ]]; then
    marcar "PING ja concluido"
    return 0
  fi
  cd "$BASE/tools/PING"
  apptainer exec --bind "$BASE" ping.sif Rscript PING_run.R \
    --fqDirectory "$BASE/fastq" \
    --resultsDirectory "$BASE/ping" \
    --fastqPattern fq.gz \
    --threads "$THREADS" >> "$BASE/logs/ping.log" 2>&1
  cd "$BASE"

  if [[ -s $BASE/ping/finalAlleleCalls.csv ]]; then
    marcar "PING ok"
    return 0
  elif [[ -s $BASE/ping/iterAlleleCalls.csv ]]; then
    marcar "PING parou antes do final; usando iterAlleleCalls.csv"
    return 0
  else
    marcar "PING FALHOU - ver logs/ping.log (comum com poucas amostras)"
    return 1
  fi
}

# ---------------------------------------------------------------- 4. verdade
etapa_verdade() {
  titulo "4. TABELA DE VERDADE"
  [[ -s $BASE/verdade_hprc.tsv ]] && { marcar "verdade ja gerada"; return 0; }
  python3 "$BASE/immuannot_para_verdade.py" "$BASE/verdade/hprc" \
    -o "$BASE/verdade_hprc.tsv" --amostras "$BASE/amostras.txt" \
    >> "$BASE/logs/verdade.log" 2>&1
  if [[ -s $BASE/verdade_hprc.tsv ]]; then
    marcar "verdade: $(( $(wc -l < "$BASE/verdade_hprc.tsv") - 1 )) pares amostra x locus"
  else
    marcar "verdade FALHOU"
    return 1
  fi
}

# ---------------------------------------------------------------- 5. kirsolve
etapa_kirsolve() {
  titulo "5. KIRSOLVE"
  source "$BASE/tools/kirsolve/.venv/bin/activate"

  local entradas=()
  [[ -d $BASE/kirmapper/genotype/cds ]] && entradas+=("$BASE/kirmapper/genotype/cds")
  [[ -d $BASE/ping ]] && entradas+=("$BASE/ping")
  if (( ${#entradas[@]} == 0 )); then
    marcar "kirsolve: nenhuma entrada disponivel"
    return 1
  fi

  python -m kirsolve inspect "${entradas[@]}" > "$BASE/logs/inspect.txt" 2>&1
  python -m kirsolve benchmark "${entradas[@]}" \
    --truth "$BASE/verdade_hprc.tsv" \
    -o "$BASE/kirsolve" >> "$BASE/logs/benchmark.log" 2>&1

  if [[ -s $BASE/kirsolve/benchmark_resumo.json ]]; then
    marcar "benchmark ok"
    return 0
  fi
  marcar "benchmark FALHOU - ver logs/benchmark.log"
  return 1
}

# ---------------------------------------------------------------- principal
main() {
  : > "$ESTADO"
  marcar "inicio | paralelo=$PARALELO threads=$THREADS"

  titulo "0. REQUISITOS"
  if ! verificar_requisitos; then
    marcar "ABORTADO: requisitos faltando"
    exit 1
  fi
  echo "  tudo presente"

  etapa_extracao   || marcar "etapa 1 com falhas, seguindo"
  etapa_kirmapper  || marcar "etapa 2 com falhas, seguindo"
  etapa_ping       || marcar "etapa 3 com falhas, seguindo"
  etapa_verdade    || marcar "etapa 4 com falhas, seguindo"
  etapa_kirsolve   || marcar "etapa 5 com falhas"

  titulo "FIM"
  marcar "concluido"
  echo
  if [[ -f $BASE/logs/amostras_falhadas.txt ]]; then
    echo "Amostras que nao baixaram (rode o script de novo para tentar):"
    sort -u "$BASE/logs/amostras_falhadas.txt"
    echo
  fi
  [[ -s $BASE/kirsolve/benchmark_resumo.json ]] && cat "$BASE/kirsolve/benchmark_resumo.json"
  echo
  echo "Leia: $ESTADO"
  echo "      $BASE/kirsolve/benchmark_por_locus.tsv"
  echo
  echo "LEMBRETE: os limiares de copy number do ncopy foram automaticos."
  echo "Revise os graficos em kirmapper/ncopy/plots/ e, se mudar algo em"
  echo "thresholds.txt, rode ncopy e genotype de novo."
}

main "$@"
