#!/usr/bin/env bash
# Step 03: quality control of the raw reads.
#
# Runs FastQC on every raw FASTQ and summarises each dataset with MultiQC:
#   results/qc/raw/fastqc/{A,B}/        FastQC reports
#   results/qc/raw/multiqc_{A,B}.html   MultiQC reports
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/03_qc_raw.sh
set -euo pipefail

cd "$(dirname "$0")/.."

THREADS="${THREADS:-8}"
OUTDIR="results/qc/raw"

for dataset in A B; do
    echo "FastQC on dataset ${dataset}"
    mkdir -p "${OUTDIR}/fastqc/${dataset}"
    fastqc --quiet --threads "${THREADS}" --outdir "${OUTDIR}/fastqc/${dataset}" \
        data/raw/"${dataset}"/*.fastq.gz

    echo "MultiQC on dataset ${dataset}"
    multiqc --force --quiet --filename "multiqc_${dataset}" --outdir "${OUTDIR}" \
        "${OUTDIR}/fastqc/${dataset}"
done
echo "Done"
