#!/usr/bin/env bash
# Step 08: differential expression with DESeq2.
#
# Dataset A: ObP1 against ObN1, run separately for each library preparation (as
# in Zheludev et al., 2024) on the featureCounts (step 06) and kallisto (step 07)
# counts, with scripts/de_deseq2.R (see its header for the decisions: obelisk
# strands added up, rRNA genes left out, host-only PCA).
#
# Dataset B: oxygen and sugar contrasts on the featureCounts counts, with
# scripts/de_conditions.R (same gene handling, design ~ oxygen + sugar). It runs
# twice: on all 10 samples, and without the two oxic replicate 3 samples
# (glucose_oxic_3, sucrose_oxic_3), a sensitivity run because they may carry a
# batch effect (22-24 % rRNA and 14-16 % obelisk reads).
#
# Outputs:
#   A: results/de/A/<prep>/deseq2_{featurecounts,kallisto}.tsv
#      results/de/A/<prep>/size_factors.tsv
#      results/de/A/<prep>/summary.tsv
#      results/de/A/<prep>/{pca,ma,volcano,lfc_concordance}.png
#   B: results/de/B/<run>/deseq2_{oxygen,sugar}.tsv, with <run> = all, no_oxic_rep3
#      results/de/B/<run>/summary.tsv
#      results/de/B/<run>/obelisk_counts.tsv
#      results/de/B/<run>/{pca,ma,volcano,obelisk_counts}.png
#
# Usage (inside the in-silico-task environment):
#   bash scripts/08_de.sh A
#   bash scripts/08_de.sh B
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 1 || ! "$1" =~ ^[AB]$ ]]; then
    echo "Usage: $0 <dataset: A or B>" >&2
    exit 1
fi

DATASET="$1"
case "${DATASET}" in
    A)
        PREPS=(totalRNA rRNAdep RNaseR)
        for prep in "${PREPS[@]}"; do
            echo "DESeq2 on dataset ${DATASET}, preparation ${prep}"
            Rscript scripts/de_deseq2.R \
                --samples "config/samples_${DATASET}.tsv" --prep "${prep}" \
                --featurecounts "results/counts/${DATASET}/featurecounts.tsv" \
                --kallisto-dir "results/quant/${DATASET}" \
                --genes data/ref/genes.tsv \
                --outdir "results/de/${DATASET}/${prep}"
        done
        ;;
    B)
        INPUTS=(
            --samples config/samples_B.tsv
            --featurecounts results/counts/B/featurecounts.tsv
            --genes data/ref/genes.tsv
        )
        echo "DESeq2 on dataset B, all samples"
        Rscript scripts/de_conditions.R "${INPUTS[@]}" --outdir results/de/B/all
        echo "DESeq2 on dataset B, without the oxic replicate 3 samples"
        Rscript scripts/de_conditions.R "${INPUTS[@]}" \
            --exclude glucose_oxic_3,sucrose_oxic_3 \
            --outdir results/de/B/no_oxic_rep3
        ;;
esac
echo "Done"
