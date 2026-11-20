#!/usr/bin/env bash
# Step 09: compare the obelisk metrics of this work with the published values.
#
# Dataset A is compared with Zheludev et al. (2024) per library preparation and
# dataset B with Maddamsetti & You (2025), per oxygen condition and per run, with
# scripts/compare_studies.py.
#
# Inputs:
#   results/obelisk/{A,B}/obelisk_metrics.tsv   our metrics (step 06)
#   resources/zheludev2024_values.tsv           published values, dataset A
#   resources/maddamsetti2025_values.tsv        published values, dataset B
#
# Outputs:
#   results/comparison/obelisk_summary.tsv      one row per group, "mean ± sd"
#   results/comparison/b_vs_maddamsetti.tsv     one row per run of dataset B
#
# Usage (inside the in-silico-task environment):
#   bash scripts/09_compare.sh
set -euo pipefail

cd "$(dirname "$0")/.."

echo "Comparison with the published values"
python scripts/compare_studies.py \
    --metrics-a results/obelisk/A/obelisk_metrics.tsv \
    --metrics-b results/obelisk/B/obelisk_metrics.tsv \
    --zheludev resources/zheludev2024_values.tsv \
    --maddamsetti resources/maddamsetti2025_values.tsv \
    --outdir results/comparison
echo "Done"
