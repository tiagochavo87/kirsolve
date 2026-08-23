#!/bin/bash
#SBATCH --job-name=kirsolve
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/kirsolve_%j.out

# Uso: sbatch scripts/slurm_kirsolve.sh dados.xlsx resultados/
# O gargalo do pipeline e memoria (expansao combinatoria dos candidatos),
# nao CPU. 16 GB cobrem coortes de milhares de amostras.

set -euo pipefail
INPUT=${1:?informe a planilha de entrada}
OUT=${2:-resultados}

module load python/3.12 2>/dev/null || true
source .venv/bin/activate 2>/dev/null || true

python -m kirsolve run "$INPUT" -o "$OUT" \
  --sheet-tool Sheet5=kir-mapper \
  --use-ipd \
  --min-samples-em 20
