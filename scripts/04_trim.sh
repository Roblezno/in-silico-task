#!/usr/bin/env bash
# Step 04: trim reads with Cutadapt and re-run quality control.
#
# Dataset-specific options, chosen from the raw FastQC reports:
#   A (SMART-Seq Total RNA Mid Input, MiSeq): remove the 3 template-switching
#     nucleotides (GGG) at the start of read 1, the Illumina adapter and
#     low-quality 3' ends (-q 20).
#   B (Illumina stranded kit, NovaSeq): remove the first base of both reads
#     (a non-biological T in 99.9% of reads), the Nextera-type adapter and
#     low-quality 3' ends with two-colour trimming (also removes poly-G tails).
# Both: discard pairs in which either read is shorter than 20 nt.
#
# Outputs:
#   data/trimmed/{A,B}/<run>_{1,2}.fastq.gz   trimmed reads
#   results/qc/trimmed/cutadapt/{A,B}/<run>.log
#   results/qc/trimmed/fastqc/{A,B}/
#   results/qc/trimmed/multiqc_{A,B}.html     FastQC + Cutadapt summary
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/04_trim.sh
set -euo pipefail

cd "$(dirname "$0")/.."

THREADS="${THREADS:-8}"
MIN_LENGTH=20
OUTDIR="results/qc/trimmed"

OPTIONS_A=(-u 3 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20)
OPTIONS_B=(-u 1 -U 1 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT --nextseq-trim=20)

for dataset in A B; do
    if [[ "${dataset}" == "A" ]]; then
        options=("${OPTIONS_A[@]}")
    else
        options=("${OPTIONS_B[@]}")
    fi
    mkdir -p "data/trimmed/${dataset}" "${OUTDIR}/cutadapt/${dataset}" "${OUTDIR}/fastqc/${dataset}"

    for run in $(tail -n +2 "config/samples_${dataset}.tsv" | cut -f1); do
        echo "Cutadapt on ${run} (dataset ${dataset})"
        cutadapt "${options[@]}" --minimum-length "${MIN_LENGTH}" --cores "${THREADS}" \
            -o "data/trimmed/${dataset}/${run}_1.fastq.gz" \
            -p "data/trimmed/${dataset}/${run}_2.fastq.gz" \
            "data/raw/${dataset}/${run}_1.fastq.gz" "data/raw/${dataset}/${run}_2.fastq.gz" \
            > "${OUTDIR}/cutadapt/${dataset}/${run}.log"
    done

    echo "FastQC on trimmed dataset ${dataset}"
    fastqc --quiet --threads "${THREADS}" --outdir "${OUTDIR}/fastqc/${dataset}" \
        data/trimmed/"${dataset}"/*.fastq.gz

    echo "MultiQC on trimmed dataset ${dataset}"
    multiqc --force --quiet --filename "multiqc_${dataset}" --outdir "${OUTDIR}" \
        "${OUTDIR}/cutadapt/${dataset}" "${OUTDIR}/fastqc/${dataset}"
done
echo "Done"
