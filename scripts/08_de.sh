#!/usr/bin/env bash
# Step 08: differential expression of ObP1 against ObN1 with DESeq2, run
# separately for each library preparation (as in Zheludev et al., 2024) on the
# featureCounts (step 06) and kallisto (step 07) counts.
#
# The analysis itself is in scripts/de_deseq2.R (see its header for the
# decisions: obelisk strands added up, rRNA genes left out, host-only PCA).
#
# Outputs:
#   results/de/<dataset>/<prep>/deseq2_{featurecounts,kallisto}.tsv
#   results/de/<dataset>/<prep>/size_factors.tsv
#   results/de/<dataset>/<prep>/summary.tsv
#   results/de/<dataset>/<prep>/{pca,ma,volcano,lfc_concordance}.png
#
# Usage (inside the in-silico-task environment):
#   bash scripts/08_de.sh A
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 1 || ! "$1" =~ ^[AB]$ ]]; then
    echo "Usage: $0 <dataset: A or B>" >&2
    exit 1
fi

DATASET="$1"
case "${DATASET}" in
    A) PREPS=(totalRNA rRNAdep RNaseR) ;;
    B) echo "The DE design of dataset B is not set yet (phase F5)" >&2; exit 1 ;;
esac

for prep in "${PREPS[@]}"; do
    echo "DESeq2 on dataset ${DATASET}, preparation ${prep}"
    Rscript scripts/de_deseq2.R \
        --samples "config/samples_${DATASET}.tsv" --prep "${prep}" \
        --featurecounts "results/counts/${DATASET}/featurecounts.tsv" \
        --kallisto-dir "results/quant/${DATASET}" \
        --genes data/ref/genes.tsv \
        --outdir "results/de/${DATASET}/${prep}"
done
echo "Done"
