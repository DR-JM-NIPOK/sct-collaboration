#!/bin/bash
# run_chains.sh - Run nested sampling chains for CAR analysis
# Usage: bash run_chains.sh --data combined --model car

set -e

# Default values
DATA="combined"
MODEL="car"
OUTPUT_DIR="/output"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data)
            DATA="$2"
            shift 2
            ;;
        --model)
            MODEL="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "========================================="
echo "SCT Cosmology Series - Chain Runner"
echo "========================================="
echo "Data: $DATA"
echo "Model: $MODEL"
echo "Output: $OUTPUT_DIR"
echo "========================================="

# Run PolyChord
pypolychord \
    --config chains/config_${MODEL}.ini \
    --data chains/data_${DATA}.txt \
    --output ${OUTPUT_DIR}/chains_${MODEL}_${DATA}

# Compute evidence
python chains/compute_evidence.py \
    --input ${OUTPUT_DIR}/chains_${MODEL}_${DATA} \
    --output ${OUTPUT_DIR}/evidence_${MODEL}_${DATA}.txt

# Generate posterior plots
python chains/plot_posteriors.py \
    --input ${OUTPUT_DIR}/chains_${MODEL}_${DATA} \
    --output ${OUTPUT_DIR}/posteriors_${MODEL}_${DATA}.png

echo "Done."