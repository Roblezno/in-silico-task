#!/usr/bin/env bash
# Step 10: benchmark of four obelisk reference variants (R0-R3) x five mapping
# configurations on one sample per dataset.
#
# - One sample per dataset is enough: the aim is to compare how references and
#   mappers treat obelisk reads (mapped fraction, junction reads, speed), not
#   to compare biological conditions, so replicates add nothing and 20
#   combinations per dataset stay quick.
# - Reference variants (same host genome, obelisk extended by N nt):
#     R0  N = 0                linear obelisk
#     R1  N = 30               the reference of the main analysis
#     R2  N = max read length - 1 (A: 77, B: 98)
#     R3  N = 1137             obelisk length: the contig is a dimer
# - Configurations: bwa mem, bowtie2 end-to-end, bowtie2 --local, HISAT2
#   (no spliced alignment) and kallisto, all with default parameters.
# - Dataset A uses SRR29668056 (ObP1-A_totalRNA) whole. Dataset B uses
#   SRR23591557 (sucrose, oxic, replicate 1) subsampled to 2 M pairs.
# - Unmapped reads stay in the BAM, so samtools flagstat reports the total
#   number of read pairs.
# - Only the mapper process is timed (GNU time), not samtools sort nor the
#   index building. References, indexes and mappings are done once: a finished
#   mapping leaves a <ref>.<tool>.done marker in results/benchmark/<dataset>/logs/
#   and is skipped on the next run. Delete the markers to re-map (the time
#   measurements come from the run that made the BAM).
# - Counting and metrics are quick and always re-run: featureCounts per
#   variant (primary alignments, multi-mapping pairs included), the obelisk
#   metrics of the main analysis (scripts/obelisk_metrics.py) and the
#   benchmark summary (scripts/benchmark_metrics.py).
#
# Outputs:
#   data/benchmark/B/<run>_sub_{1,2}.fastq.gz           subsampled reads (B only)
#   data/ref/benchmark/<dataset>/<ref>/                 reference files of the
#                                                       variant and its indexes
#                                                       (bwa/, bowtie2/, hisat2/,
#                                                       kallisto/)
#   data/benchmark/<dataset>/bam/<ref>.<tool>.bam(.bai) coordinate-sorted alignments
#   results/benchmark/<dataset>/input_pairs.tsv         read pairs given to the mappers
#   results/benchmark/<dataset>/logs/<ref>.<tool>.log         mapper stderr
#   results/benchmark/<dataset>/logs/<ref>.<tool>.time.txt    GNU time -v report
#   results/benchmark/<dataset>/logs/<ref>.<tool>.flagstat    (aligners only)
#   results/benchmark/<dataset>/logs/<ref>.<tool>.done        mapping finished marker
#   results/benchmark/<dataset>/kallisto/<ref>/         abundance.tsv, run_info.json
#   results/benchmark/<dataset>/counts/<ref>.featurecounts.tsv(.summary,.log)
#                                                       counts of the four aligners
#   results/benchmark/<dataset>/metrics/<ref>.combinations.tsv     sample sheet of
#                                                       the aligners (run = <ref>.<tool>)
#   results/benchmark/<dataset>/metrics/<ref>.obelisk_metrics.tsv  obelisk metrics
#   results/benchmark/<dataset>/benchmark_summary.tsv   one row per variant x tool
#   results/benchmark/<dataset>/host_spearman_R1.tsv    Spearman correlation of the
#                                                       host gene counts between tools
#   <tool> is bwa_mem, bowtie2_e2e, bowtie2_local, hisat2 or kallisto.
#
# Usage (inside the in-silico-task environment):
#   THREADS=8 bash scripts/10_benchmark.sh A
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 1 || ! "$1" =~ ^[AB]$ ]]; then
    echo "Usage: $0 <dataset: A or B>" >&2
    exit 1
fi

DATASET="$1"
THREADS="${THREADS:-8}"

# Sample, reads to keep (0 = whole sample), R2 extension and strandedness.
case "${DATASET}" in
    A)
        RUN="SRR29668056"       # ObP1-A_totalRNA, ~0.70 M pairs
        SUBSAMPLE=0
        R2_EXTENSION=77         # trimmed A reads are at most 78 nt
        STRAND_OPTION="--fr-stranded"   # read 1 has the RNA strand (featureCounts -s 1)
        STRAND=1
        ;;
    B)
        RUN="SRR23591557"       # sucrose, oxic, replicate 1
        SUBSAMPLE=2000000
        R2_EXTENSION=98         # trimmed B reads are at most 99 nt
        STRAND_OPTION="--rf-stranded"   # dUTP, read 1 is antisense (featureCounts -s 2)
        STRAND=2
        ;;
esac

GENBANK="data/ref/CP000387.1.gb"
OBELISK_FASTA="resources/obelisk_000003.fasta"
REFROOT="data/ref/benchmark/${DATASET}"
BAMDIR="data/benchmark/${DATASET}/bam"
OUTDIR="results/benchmark/${DATASET}"
LOGDIR="${OUTDIR}/logs"

# Variant name -> nt of the obelisk start appended to its end.
VARIANTS=(R0 R1 R2 R3)
declare -A EXTENSION=([R0]=0 [R1]=30 [R2]="${R2_EXTENSION}" [R3]=1137)
INDEX_TOOLS=(bwa bowtie2 hisat2 kallisto)
ALIGNER_TOOLS=(bwa_mem bowtie2_e2e bowtie2_local hisat2)
MAP_TOOLS=("${ALIGNER_TOOLS[@]}" kallisto)

# --- Reads ---------------------------------------------------------------

TRIMMED="data/trimmed/${DATASET}/${RUN}"
if [[ "${SUBSAMPLE}" -eq 0 ]]; then
    READ1="${TRIMMED}_1.fastq.gz"
    READ2="${TRIMMED}_2.fastq.gz"
else
    READDIR="data/benchmark/${DATASET}"
    READ1="${READDIR}/${RUN}_sub_1.fastq.gz"
    READ2="${READDIR}/${RUN}_sub_2.fastq.gz"
    # The marker file is written only after both files are complete.
    if [[ ! -f "${READDIR}/${RUN}_sub.done" ]]; then
        echo "Subsampling ${RUN} to ${SUBSAMPLE} pairs"
        mkdir -p "${READDIR}"
        # The same seed on both mates selects the same records, so the pairs
        # stay paired.
        seqtk sample -s 42 "${TRIMMED}_1.fastq.gz" "${SUBSAMPLE}" | gzip > "${READ1}"
        seqtk sample -s 42 "${TRIMMED}_2.fastq.gz" "${SUBSAMPLE}" | gzip > "${READ2}"
        touch "${READDIR}/${RUN}_sub.done"
    fi
fi

mkdir -p "${BAMDIR}" "${LOGDIR}" "${OUTDIR}/kallisto"
lines=$(gzip -dc "${READ1}" | wc -l)
PAIRS=$((lines / 4))
printf 'run\tpairs\n%s\t%d\n' "${RUN}" "${PAIRS}" > "${OUTDIR}/input_pairs.tsv"

# --- References and indexes ---------------------------------------------

# Build the index of one tool for one variant unless its marker file exists.
# The marker is written only after a successful build, so an interrupted build
# is not mistaken for a complete index.
build_index() {
    local ref="$1" tool="$2"
    local dir="${REFROOT}/${ref}"
    local prefix="${dir}/${tool}/index"
    if [[ -f "${prefix}.done" ]]; then
        return
    fi
    echo "Building ${tool} index for ${ref}"
    mkdir -p "${dir}/${tool}"
    case "${tool}" in
        bwa) bwa index -p "${prefix}" "${dir}/sk36_obelisk.fa" ;;
        bowtie2) bowtie2-build --threads "${THREADS}" "${dir}/sk36_obelisk.fa" "${prefix}" ;;
        hisat2) hisat2-build -p "${THREADS}" "${dir}/sk36_obelisk.fa" "${prefix}" ;;
        kallisto)
            kallisto index -t "${THREADS}" -T "${dir}/kallisto/tmp" \
                -i "${prefix}.idx" "${dir}/transcripts.fa"
            ;;
    esac > "${prefix}.build.log" 2>&1
    touch "${prefix}.done"
}

for ref in "${VARIANTS[@]}"; do
    dir="${REFROOT}/${ref}"
    if [[ ! -f "${dir}/reference.done" ]]; then
        echo "Building reference ${ref} (extension ${EXTENSION[${ref}]} nt)"
        mkdir -p "${dir}"
        python scripts/build_references.py \
            --genbank "${GENBANK}" \
            --obelisk "${OBELISK_FASTA}" \
            --extension "${EXTENSION[${ref}]}" \
            --outdir "${dir}"
        samtools faidx "${dir}/sk36_obelisk.fa"
        touch "${dir}/reference.done"
    fi
    for tool in "${INDEX_TOOLS[@]}"; do
        build_index "${ref}" "${tool}"
    done
done

# --- Mapping -------------------------------------------------------------

# Run a command under GNU time, writing the report to the file in $1. `time`
# is a bash keyword, so `command` is needed to reach the GNU binary.
timed() {
    local report="$1"
    shift
    command time -v -o "${report}" "$@"
}

# Align with one of the four aligners, then sort, index and count the BAM.
# The mapper stderr (alignment summary) goes to its log file.
map_aligner() {
    local ref="$1" tool="$2"
    local dir="${REFROOT}/${ref}"
    local prefix="${LOGDIR}/${ref}.${tool}"
    local bam="${BAMDIR}/${ref}.${tool}.bam"
    local -a cmd
    case "${tool}" in
        bwa_mem)
            cmd=(bwa mem -t "${THREADS}" "${dir}/bwa/index" "${READ1}" "${READ2}")
            ;;
        bowtie2_e2e)
            cmd=(bowtie2 -p "${THREADS}" -x "${dir}/bowtie2/index" -1 "${READ1}" -2 "${READ2}")
            ;;
        bowtie2_local)
            cmd=(bowtie2 --local -p "${THREADS}" -x "${dir}/bowtie2/index" -1 "${READ1}" -2 "${READ2}")
            ;;
        hisat2)
            cmd=(hisat2 -p "${THREADS}" --no-spliced-alignment -x "${dir}/hisat2/index"
                -1 "${READ1}" -2 "${READ2}")
            ;;
    esac
    timed "${prefix}.time.txt" "${cmd[@]}" 2> "${prefix}.log" \
        | samtools sort -@ "${THREADS}" -o "${bam}" -
    samtools index "${bam}"
    samtools flagstat -@ "${THREADS}" "${bam}" > "${prefix}.flagstat"
}

# Pseudoalign with kallisto (stranded, plaintext output).
map_kallisto() {
    local ref="$1"
    local prefix="${LOGDIR}/${ref}.kallisto"
    mkdir -p "${OUTDIR}/kallisto/${ref}"
    timed "${prefix}.time.txt" \
        kallisto quant -t "${THREADS}" "${STRAND_OPTION}" --plaintext \
        -i "${REFROOT}/${ref}/kallisto/index.idx" -o "${OUTDIR}/kallisto/${ref}" \
        "${READ1}" "${READ2}" 2> "${prefix}.log"
}

# The marker is written only after the whole combination is finished (BAM,
# index and flagstat, or the kallisto output), so an interrupted mapping is
# redone.
for ref in "${VARIANTS[@]}"; do
    for tool in "${MAP_TOOLS[@]}"; do
        marker="${LOGDIR}/${ref}.${tool}.done"
        if [[ -f "${marker}" ]]; then
            echo "${tool} on ${ref} already done, skipping"
            continue
        fi
        echo "${tool} on ${ref} (dataset ${DATASET}, ${RUN})"
        if [[ "${tool}" == "kallisto" ]]; then
            map_kallisto "${ref}"
        else
            map_aligner "${ref}" "${tool}"
        fi
        touch "${marker}"
    done
done

# --- Counting and metrics ------------------------------------------------

COUNTDIR="${OUTDIR}/counts"
METRICDIR="${OUTDIR}/metrics"
mkdir -p "${COUNTDIR}" "${METRICDIR}"

# Per variant: count the four aligners together, then compute the obelisk
# metrics of the main analysis. Always re-run (quick).
for ref in "${VARIANTS[@]}"; do
    dir="${REFROOT}/${ref}"
    sheet="${METRICDIR}/${ref}.combinations.tsv"
    counts="${COUNTDIR}/${ref}.featurecounts.tsv"

    # obelisk_metrics.py reads <bamdir>/<run>.bam and featureCounts names its
    # columns by the BAM stem, so the combination name <ref>.<tool> plays the
    # role of the run. The script drops read_count.
    {
        printf 'run\tref\ttool\tread_count\n'
        for tool in "${ALIGNER_TOOLS[@]}"; do
            printf '%s.%s\t%s\t%s\t%d\n' "${ref}" "${tool}" "${ref}" "${tool}" "${PAIRS}"
        done
    } > "${sheet}"

    bams=()
    for tool in "${ALIGNER_TOOLS[@]}"; do
        bams+=("${BAMDIR}/${ref}.${tool}.bam")
    done

    # With --primary every primary alignment is counted, also of multi-mapping
    # pairs (R2 and R3 repeat obelisk sequence), so each pair counts once at
    # its primary alignment, the same rule as obelisk_metrics.py. bwa and
    # bowtie2 set no NH tag, so this is the only rule that treats all aligners
    # alike.
    echo "featureCounts on ${ref} (dataset ${DATASET}, -s ${STRAND})"
    featureCounts -p --countReadPairs -s "${STRAND}" -M --primary -T "${THREADS}" \
        -F SAF -a "${dir}/annotation.saf" \
        -o "${counts}" "${bams[@]}" \
        > "${COUNTDIR}/${ref}.featurecounts.log" 2>&1

    echo "Obelisk metrics on ${ref} (dataset ${DATASET})"
    python scripts/obelisk_metrics.py \
        --samples "${sheet}" --bamdir "${BAMDIR}" --genes "${dir}/genes.tsv" \
        --counts "${counts}" --strand "${STRAND}" --threads "${THREADS}" \
        --output "${METRICDIR}/${ref}.obelisk_metrics.tsv"
done

echo "Benchmark summary on dataset ${DATASET}"
python scripts/benchmark_metrics.py \
    --bamdir "${BAMDIR}" --refdir "${REFROOT}" --resultdir "${OUTDIR}" \
    --obelisk "${OBELISK_FASTA}" --threads "${THREADS}"
echo "Done"
