#!/usr/bin/env bash
# Step 01: prepare the reference files for the host genome and Obelisk-S.s.
#
# Downloads the S. sanguinis SK36 genome (GenBank CP000387.1) from NCBI and
# builds, in data/ref/, the genome + obelisk FASTA (with its samtools index),
# the featureCounts SAF, the RSeQC BED12, the kallisto transcripts and a gene
# table. See scripts/build_references.py for details on each file.
#
# Usage (inside the in-silico-task environment):
#   bash scripts/01_prepare_references.sh
set -euo pipefail

cd "$(dirname "$0")/.."

ACCESSION="CP000387.1"
OBELISK_FASTA="resources/obelisk_000003.fasta"
# nt of the obelisk start appended to its end to cover the circle junction,
# as in Zheludev et al. (2024).
OBELISK_EXTENSION=30
OUTDIR="data/ref"
GENBANK="${OUTDIR}/${ACCESSION}.gb"
EFETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

mkdir -p "${OUTDIR}"

if [[ ! -s "${GENBANK}" ]]; then
    echo "Downloading ${ACCESSION} from NCBI"
    curl -fsSL "${EFETCH_URL}?db=nuccore&id=${ACCESSION}&rettype=gbwithparts&retmode=text" \
        -o "${GENBANK}"
fi

echo "Building reference files in ${OUTDIR}"
python scripts/build_references.py \
    --genbank "${GENBANK}" \
    --obelisk "${OBELISK_FASTA}" \
    --extension "${OBELISK_EXTENSION}" \
    --outdir "${OUTDIR}"

samtools faidx "${OUTDIR}/sk36_obelisk.fa"
echo "Done"
