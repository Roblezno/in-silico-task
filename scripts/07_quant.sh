#!/usr/bin/env bash
# Step 07: quantify transcripts with kallisto, as the second route for the
# differential expression analysis (step 08).
#
# - The index is built from data/ref/transcripts.fa: every host gene plus the
#   obelisk sense and antisense strands as separate targets (extended 30 nt).
#   The two obelisk targets share all their k-mers, so only stranded
#   quantification can tell them apart.
# - Strandedness, set per dataset from RSeQC infer_experiment (step 05):
#     A: read 1 has the strand of the RNA, so --fr-stranded (featureCounts -s 1).
#     B: not set yet (phase F5).
# - No bootstraps: DESeq2 does not use them.
# - Only obelisk-positive and obelisk-negative runs are quantified; media
#   controls have too few reads and are not part of the DE analysis.
# - --plaintext writes abundance.tsv instead of HDF5, so tximport does not need
#   rhdf5.
#
# Outputs:
#   data/ref/kallisto/transcripts.idx                 kallisto index (built once)
#   results/quant/<dataset>/<run>/abundance.tsv       estimated counts and TPM
#   results/quant/<dataset>/<run>/run_info.json       pseudoalignment summary
#   results/quant/<dataset>/<run>/kallisto.log
#   results/quant/multiqc_quant_<dataset>.html
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/07_quant.sh A
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 1 || ! "$1" =~ ^[AB]$ ]]; then
    echo "Usage: $0 <dataset: A or B>" >&2
    exit 1
fi

DATASET="$1"
THREADS="${THREADS:-8}"
case "${DATASET}" in
    A) STRAND_OPTION="--fr-stranded" ;;
    B) echo "Strandedness of dataset B is not set yet (phase F5)" >&2; exit 1 ;;
esac
SAMPLES="config/samples_${DATASET}.tsv"
TRANSCRIPTS="data/ref/transcripts.fa"
INDEXDIR="data/ref/kallisto"
INDEX="${INDEXDIR}/transcripts.idx"
OUTDIR="results/quant/${DATASET}"

# The marker file is written only after a successful build, so an interrupted
# build is not mistaken for a complete index.
if [[ ! -f "${INDEX}.done" ]]; then
    echo "Building kallisto index"
    mkdir -p "${INDEXDIR}"
    kallisto index -t "${THREADS}" -T "${INDEXDIR}/tmp" -i "${INDEX}" "${TRANSCRIPTS}" \
        > "${INDEX}.build.log" 2>&1
    touch "${INDEX}.done"
fi

mkdir -p "${OUTDIR}"

# Runs of the two strains only (column 3 of the sample sheet).
runs=$(awk -F'\t' 'NR > 1 && ($3 == "ObP1" || $3 == "ObN1") {print $1}' "${SAMPLES}")

for run in ${runs}; do
    echo "kallisto on ${run} (dataset ${DATASET}, ${STRAND_OPTION})"
    mkdir -p "${OUTDIR}/${run}"
    kallisto quant -t "${THREADS}" ${STRAND_OPTION} --plaintext \
        -i "${INDEX}" -o "${OUTDIR}/${run}" \
        "data/trimmed/${DATASET}/${run}_1.fastq.gz" \
        "data/trimmed/${DATASET}/${run}_2.fastq.gz" \
        > "${OUTDIR}/${run}/kallisto.log" 2>&1
done

echo "MultiQC on dataset ${DATASET}"
multiqc --force --quiet --filename "multiqc_quant_${DATASET}" --outdir results/quant "${OUTDIR}"
echo "Done"
