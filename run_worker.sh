#!/usr/bin/env bash
# run_worker.sh — single-sequence launcher (SLURM-friendly)

#SBATCH --job-name=gbnl_one
#SBATCH --time=00:15:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
# NOTE: make sure 'logs' exists before submission
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

set -euo pipefail

# ---------------- inputs readers can edit ----------------
SEQ_ID="N-China-F"
SEQ="AGGGAACTTCTCCTGCTAGAAT" 

# Output location (relative to the repository)
OUT_ROOT="./outputs"

# Script path 
MAIN_PY="./main.py"           

# VR and filtration settings
VR_DIM=1
MIN_E=0
MAX_E=10

# Macaulay2 backend
USE_CONTAINER=true
M2_IMG="${HOME}/m2-1.24.05.sif"

# Python interpreter
PYTHON_BIN="${PYTHON_BIN:-python}"

# ----------------------------------------------------------
# IMPORTANT: create 'logs' *before* submitting, e.g.:
#   mkdir -p logs
# or use absolute paths in #SBATCH --output/--error.

mkdir -p "$OUT_ROOT"

# Normalize sequence for logs
SEQ_UPPER=$(echo -n "$SEQ" | tr '[:lower:]' '[:upper:]' | sed 's/U/T/g')

echo "GBNL single-sequence run"
echo "  seq id     : $SEQ_ID"
echo "  length     : ${#SEQ_UPPER}"
echo "  out root   : $OUT_ROOT"
echo "  vr dim     : $VR_DIM"
echo "  eps range  : [$MIN_E, $MAX_E]"
echo "  container  : ${USE_CONTAINER}"
echo

run_one_nuc () {
  local nuc="$1"
  echo ">>> nucleotide = ${nuc}"
  local args=(
    --seq-id "$SEQ_ID"
    --seq-ord "1"
    --sequence "$SEQ_UPPER"
    --nuc "$nuc"
    --out-root "$OUT_ROOT"
    --vr-dim "$VR_DIM"
    --min-e "$MIN_E"
    --max-e "$MAX_E"
  )
  if [[ "$USE_CONTAINER" == "true" ]]; then
    args+=( --use-container --m2-image "$M2_IMG" )
  fi

  # Use all CPUs allocated to the task (or fall back to 1 when running locally)
  local C=${SLURM_CPUS_PER_TASK:-1}
  if command -v srun >/dev/null 2>&1; then
    srun -c "$C" "$PYTHON_BIN" "$MAIN_PY" "${args[@]}"
  else
    "$PYTHON_BIN" "$MAIN_PY" "${args[@]}"
  fi
}

# Run all four nucleotides
for N in A C G T; do
  run_one_nuc "$N"
done

echo
echo "Done. Inspect outputs under: ${OUT_ROOT}/01_${SEQ_ID} or ${OUT_ROOT}/${SEQ_ID} (depending on seq-ord usage)."
