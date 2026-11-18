#!/usr/bin/env bash
# Step 05: map trimmed reads to the SK36 + obelisk reference with HISAT2 and
# infer the library strandedness with RSeQC.
#
# - Bacterial genes have no introns, so spliced alignment is disabled.
# - Unmapped reads are kept in the BAM, so samtools flagstat also reports the
#   total number of read pairs.
# - The strandedness is inferred on host genes only (data/ref/host_genes.bed):
#   the obelisk is very abundant and has its own strand bias. Media controls
#   (dataset A) are skipped because they have too few reads.
#
# Outputs:
#   data/ref/hisat2/sk36_obelisk.*.ht2               HISAT2 index (built once)
#   data/bam/<dataset>/<run>.bam(.bai)               coordinate-sorted alignments
#   results/mapping/<dataset>/<run>.hisat2.log       HISAT2 alignment summary
#   results/mapping/<dataset>/<run>.flagstat
#   results/mapping/<dataset>/<run>.idxstats
#   results/mapping/<dataset>/<run>.infer_experiment.txt
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/05_map.sh A
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 1 || ! "$1" =~ ^[AB]$ ]]; then
    echo "Usage: $0 <dataset: A or B>" >&2
    exit 1
fi

DATASET="$1"
THREADS="${THREADS:-8}"
SAMPLES="config/samples_${DATASET}.tsv"
REFERENCE="data/ref/sk36_obelisk.fa"
INDEX="data/ref/hisat2/sk36_obelisk"
HOST_BED="data/ref/host_genes.bed"
BAMDIR="data/bam/${DATASET}"
OUTDIR="results/mapping/${DATASET}"

# The marker file is written only after a successful build, so an interrupted
# build is not mistaken for a complete index.
if [[ ! -f "${INDEX}.done" ]]; then
    echo "Building HISAT2 index"
    mkdir -p "$(dirname "${INDEX}")"
    hisat2-build -p "${THREADS}" "${REFERENCE}" "${INDEX}" > "${INDEX}.build.log"
    touch "${INDEX}.done"
fi

mkdir -p "${BAMDIR}" "${OUTDIR}"

for run in $(tail -n +2 "${SAMPLES}" | cut -f1); do
    echo "HISAT2 on ${run} (dataset ${DATASET})"
    hisat2 -p "${THREADS}" --no-spliced-alignment -x "${INDEX}" \
        -1 "data/trimmed/${DATASET}/${run}_1.fastq.gz" \
        -2 "data/trimmed/${DATASET}/${run}_2.fastq.gz" \
        2> "${OUTDIR}/${run}.hisat2.log" \
        | samtools sort -@ "${THREADS}" -o "${BAMDIR}/${run}.bam" -
    samtools index "${BAMDIR}/${run}.bam"
    samtools flagstat -@ "${THREADS}" "${BAMDIR}/${run}.bam" > "${OUTDIR}/${run}.flagstat"
    samtools idxstats "${BAMDIR}/${run}.bam" > "${OUTDIR}/${run}.idxstats"

    if grep -qP "^${run}\t.*\tmedia_control\t" "${SAMPLES}"; then
        continue
    fi
    infer_experiment.py -r "${HOST_BED}" -i "${BAMDIR}/${run}.bam" \
        > "${OUTDIR}/${run}.infer_experiment.txt"
done
echo "Done"
