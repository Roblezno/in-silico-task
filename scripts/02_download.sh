#!/usr/bin/env bash
# Step 02: fetch run metadata, build the sample sheets and download the reads.
#
# 1. Saves a snapshot of the ENA run metadata of each BioProject in config/
#    (versioned; delete a file to refresh it).
# 2. Builds config/samples_A.tsv and config/samples_B.tsv.
# 3. Downloads every run from NCBI SRA (prefetch), checks its integrity
#    (vdb-validate) and writes gzipped paired FASTQ files to data/raw/{A,B}/.
#    ENA serves the same reads, but its download speed was too low to use.
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/02_download.sh
set -euo pipefail

cd "$(dirname "$0")/.."

THREADS="${THREADS:-8}"
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS="run_accession,sample_title,library_name,library_strategy,library_layout,instrument_model,read_count,base_count"
PROJECT_A="PRJNA1129866"
PROJECT_B="PRJNA937727"
SRA_DIR="data/sra"

for project in "${PROJECT_A}" "${PROJECT_B}"; do
    metadata="config/ena_${project}.tsv"
    if [[ ! -s "${metadata}" ]]; then
        echo "Fetching ENA metadata for ${project}"
        curl -fsSL "${ENA_API}?accession=${project}&result=read_run&fields=${ENA_FIELDS}&format=tsv" \
            -o "${metadata}"
    fi
done

echo "Building sample sheets"
python scripts/make_samplesheets.py \
    --ena-a "config/ena_${PROJECT_A}.tsv" \
    --ena-b "config/ena_${PROJECT_B}.tsv" \
    --outdir config

# A run is complete when both gzipped mates exist and no uncompressed file is left
# (gzip removes the .fastq only after it finishes).
run_is_done() {
    local prefix="$1"
    [[ -s "${prefix}_1.fastq.gz" && -s "${prefix}_2.fastq.gz" \
        && ! -e "${prefix}_1.fastq" && ! -e "${prefix}_2.fastq" ]]
}

for dataset in A B; do
    outdir="data/raw/${dataset}"
    mkdir -p "${outdir}"
    for run in $(tail -n +2 "config/samples_${dataset}.tsv" | cut -f1); do
        if run_is_done "${outdir}/${run}"; then
            continue
        fi
        echo "Downloading ${run} (dataset ${dataset})"
        prefetch --max-size u --output-directory "${SRA_DIR}" "${run}"
        vdb-validate "${SRA_DIR}/${run}/${run}.sra"
        fasterq-dump --split-3 --force --threads "${THREADS}" --temp "${SRA_DIR}" \
            --outdir "${outdir}" "${SRA_DIR}/${run}/${run}.sra"
        gzip -f "${outdir}/${run}_1.fastq" &
        gzip_pid=$!
        gzip -f "${outdir}/${run}_2.fastq"
        wait "${gzip_pid}"
    done
done
echo "Done"
