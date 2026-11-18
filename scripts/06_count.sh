#!/usr/bin/env bash
# Step 06: count read pairs per gene with featureCounts, compute the obelisk
# metrics and plot the obelisk coverage by strand.
#
# - Strandedness (-s), set per dataset from RSeQC infer_experiment (step 05):
#     A: 0.965-0.994 of the reads are explained by "1++,1--,2+-,2-+", so read 1
#        has the strand of the RNA: -s 1.
#     B: not set yet (phase F5).
# - featureCounts default options count only uniquely mapped pairs, so pairs
#   from the four identical rRNA operons are left unassigned. The rRNA fraction
#   is computed by scripts/obelisk_metrics.py instead.
# - The obelisk is annotated on both strands, so stranded counting gives one
#   count for the sense strand and one for the antisense strand.
#
# Outputs:
#   results/counts/<dataset>/featurecounts.tsv(.summary)   counts per gene
#   results/obelisk/<dataset>/obelisk_metrics.tsv          one row per run
#   results/obelisk/<dataset>/obelisk_coverage.{png,tsv}   coverage along the obelisk by strand
#   results/mapping/multiqc_<dataset>.html                 HISAT2, samtools, RSeQC and featureCounts
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/06_count.sh A
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 1 || ! "$1" =~ ^[AB]$ ]]; then
    echo "Usage: $0 <dataset: A or B>" >&2
    exit 1
fi

DATASET="$1"
THREADS="${THREADS:-8}"
case "${DATASET}" in
    A)
        STRAND=1
        # Coverage figure: obelisk-positive runs, one panel per preparation.
        COVERAGE_OPTIONS=(--keep strain=ObP1 --group-by prep --groups totalRNA rRNAdep RNaseR)
        ;;
    B) echo "Strandedness of dataset B is not set yet (phase F5)" >&2; exit 1 ;;
esac
SAMPLES="config/samples_${DATASET}.tsv"
BAMDIR="data/bam/${DATASET}"
COUNTDIR="results/counts/${DATASET}"
OBELISKDIR="results/obelisk/${DATASET}"

mkdir -p "${COUNTDIR}" "${OBELISKDIR}"

bams=()
for run in $(tail -n +2 "${SAMPLES}" | cut -f1); do
    bams+=("${BAMDIR}/${run}.bam")
done

echo "featureCounts on dataset ${DATASET} (${#bams[@]} BAM files, -s ${STRAND})"
featureCounts -p --countReadPairs -s "${STRAND}" -T "${THREADS}" \
    -F SAF -a data/ref/annotation.saf \
    -o "${COUNTDIR}/featurecounts.tsv" "${bams[@]}" \
    > "${COUNTDIR}/featurecounts.log" 2>&1

echo "Obelisk metrics on dataset ${DATASET}"
python scripts/obelisk_metrics.py \
    --samples "${SAMPLES}" --bamdir "${BAMDIR}" --strand "${STRAND}" \
    --genes data/ref/genes.tsv --counts "${COUNTDIR}/featurecounts.tsv" \
    --threads "${THREADS}" --output "${OBELISKDIR}/obelisk_metrics.tsv"

echo "Obelisk coverage on dataset ${DATASET}"
python scripts/plot_obelisk_coverage.py \
    --samples "${SAMPLES}" --bamdir "${BAMDIR}" --strand "${STRAND}" \
    --obelisk resources/obelisk_000003.fasta --metrics "${OBELISKDIR}/obelisk_metrics.tsv" \
    "${COVERAGE_OPTIONS[@]}" \
    --threads "${THREADS}" --output-prefix "${OBELISKDIR}/obelisk_coverage"

echo "MultiQC on dataset ${DATASET}"
multiqc --force --quiet --filename "multiqc_${DATASET}" --outdir results/mapping \
    "results/mapping/${DATASET}" "${COUNTDIR}"
echo "Done"
