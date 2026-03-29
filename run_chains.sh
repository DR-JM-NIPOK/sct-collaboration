#!/usr/bin/env bash
# chains/run_chains.sh
# =====================
# Master script for nested sampling analysis of the CAR model.
# Runs PolyChord on combined or individual datasets and computes
# Bayesian evidence for model comparison vs ΛCDM.
#
# Usage:
#   bash chains/run_chains.sh                          # full combined run
#   bash chains/run_chains.sh --data desi --model car  # single dataset
#   bash chains/run_chains.sh --data combined --model lcdm  # ΛCDM comparison
#
# Author : DR JM NIPOK | License: GPL-3.0

set -euo pipefail

# ─── Defaults ─────────────────────────────────────────────────────────────────
DATA="combined"
MODEL="car"
OUTPUT_DIR="${OUTPUT_DIR:-./output}"
N_LIVE=500
N_REPEATS=10
PRECISION=0.01
NTHREADS="${NTHREADS:-4}"

# ─── Parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --data)      DATA="$2";       shift 2 ;;
    --model)     MODEL="$2";      shift 2 ;;
    --output)    OUTPUT_DIR="$2"; shift 2 ;;
    --nlive)     N_LIVE="$2";     shift 2 ;;
    --threads)   NTHREADS="$2";   shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

echo "════════════════════════════════════════════════"
echo "  SCT Cosmology — PolyChord Nested Sampling"
echo "  Model  : ${MODEL}"
echo "  Data   : ${DATA}"
echo "  Output : ${OUTPUT_DIR}"
echo "  n_live : ${N_LIVE}"
echo "════════════════════════════════════════════════"

mkdir -p "${OUTPUT_DIR}/chains_${MODEL}_${DATA}"
mkdir -p "${OUTPUT_DIR}/figures"

# ─── Step 1: Run PolyChord nested sampling ────────────────────────────────────
echo ""
echo "Step 1/4: Running PolyChord..."
python chains/run_polychord.py \
  --model   "${MODEL}" \
  --data    "${DATA}" \
  --output  "${OUTPUT_DIR}/chains_${MODEL}_${DATA}" \
  --nlive   "${N_LIVE}" \
  --repeats "${N_REPEATS}" \
  --precision "${PRECISION}" \
  --threads "${NTHREADS}"

# ─── Step 2: Extract Bayesian evidence ────────────────────────────────────────
echo ""
echo "Step 2/4: Computing Bayesian evidence..."
python chains/compute_evidence.py \
  --input  "${OUTPUT_DIR}/chains_${MODEL}_${DATA}" \
  --output "${OUTPUT_DIR}/evidence_${MODEL}_${DATA}.txt" \
  --model  "${MODEL}"

# ─── Step 3: Compare with ΛCDM (if running CAR) ──────────────────────────────
if [[ "${MODEL}" == "car" ]]; then
  echo ""
  echo "Step 3/4: Computing model comparison Δln B..."
  if [[ -f "${OUTPUT_DIR}/evidence_lcdm_${DATA}.txt" ]]; then
    python chains/compute_evidence.py \
      --compare-car  "${OUTPUT_DIR}/evidence_car_${DATA}.txt" \
      --compare-lcdm "${OUTPUT_DIR}/evidence_lcdm_${DATA}.txt" \
      --output       "${OUTPUT_DIR}/bayes_factor_${DATA}.txt"
  else
    echo "  (ΛCDM chains not found; skipping comparison. Run with --model lcdm first.)"
  fi
else
  echo "Step 3/4: Skipped (not a CAR run)."
fi

# ─── Step 4: Plot posteriors ──────────────────────────────────────────────────
echo ""
echo "Step 4/4: Plotting posteriors..."
python chains/plot_posteriors.py \
  --input  "${OUTPUT_DIR}/chains_${MODEL}_${DATA}" \
  --output "${OUTPUT_DIR}/figures/posteriors_${MODEL}_${DATA}.pdf" \
  --model  "${MODEL}"

echo ""
echo "════════════════════════════════════════════════"
echo "  Analysis complete. Results in ${OUTPUT_DIR}/"
echo "  Key files:"
echo "    evidence_${MODEL}_${DATA}.txt"
echo "    figures/posteriors_${MODEL}_${DATA}.pdf"
echo "════════════════════════════════════════════════"
