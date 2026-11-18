# in-silico-task

RNA-seq reanalysis of Obelisk-S.s in *Streptococcus sanguinis* SK36 and a benchmark of read
mappers for its quantification. Final project for the course *In silico studies in
Biomedicine* (MSc in Bioinformatics).

## Data

- **PRJNA1129866**: RNA-seq from Zheludev et al. (2024), substrains ObP1 and ObN1.
- **PRJNA937727**: RNA-seq of SK36 in monoculture (oxygen × sugar design).

Data are not versioned: the scripts download them into `data/`.

## Environment

```bash
mamba env create -f environment.yml --strict-channel-priority
mamba activate in-silico-task
```

`environment.lock.yml` records the exact package versions used.

## Layout

```
config/      sample sheets and parameters
resources/   small input files that cannot be downloaded by the scripts
scripts/     numbered scripts are analysis steps, in execution order;
             unnumbered scripts are helpers called by them
report/      Quarto report
data/        downloaded data and references (not versioned)
results/     intermediate and final results (not versioned)
```

## Usage

Run every step from the repository root inside the `in-silico-task` environment.

### 01. References

```bash
bash scripts/01_prepare_references.sh
```

Downloads the *S. sanguinis* SK36 genome (GenBank CP000387.1) from NCBI and writes to
`data/ref/`:

| File | Content | Used by |
|---|---|---|
| `sk36_obelisk.fa` (+ `.fai`) | host genome + Obelisk-S.s extended by 30 nt | mappers, IGV |
| `annotation.saf` | host genes + obelisk (sense and antisense) | featureCounts |
| `host_genes.bed` | BED12, host genes only | RSeQC `infer_experiment.py` |
| `transcripts.fa` | one sequence per gene on its own strand + obelisk (sense and antisense) | kallisto |
| `genes.tsv` | gene ID, biotype, name, product, coordinates, length | metrics, report |

**Obelisk-S.s sequence.** `resources/obelisk_000003.fasta` holds Obelisk_000003 (1137 nt),
extracted once from the supplementary table `mmc2.xlsx` of Zheludev et al. (2024) with
`scripts/extract_obelisk.py` (row `Obelisk_000003_000001_000001`). Its Oblin-1 ORF lies on
the + strand (nt 51–743), which defines the obelisk *sense* strand. The obelisk is circular:
its first 30 nt are appended to its end so that reads crossing the circle junction can align,
as in the original article. The obelisk is annotated twice, as `Obelisk_000003` (+ strand,
sense) and `Obelisk_000003_antisense` (− strand), so that stranded counting gives one count
per strand.

### 02. Sample sheets and reads

```bash
THREADS=8 bash scripts/02_download.sh
```

Saves the ENA run metadata of both BioProjects (`config/ena_PRJNA*.tsv`), builds the sample
sheets and downloads the reads:

- `config/samples_A.tsv`: the 27 RNA-seq runs of PRJNA1129866 (ObP1, ObN1 and media controls ×
  replicates A–C × total RNA, rRNA-depleted and RNaseR-treated). DNA-seq and nanopore runs are left out.
- `config/samples_B.tsv`: the 10 *S. sanguinis* monoculture runs of PRJNA937727 (glucose or
  sucrose × oxic or anoxic). Co-culture runs are left out.

Reads are downloaded from NCBI SRA (`prefetch`), checked with `vdb-validate` and converted to
gzipped paired FASTQ files in `data/raw/{A,B}/` (`fasterq-dump`). ENA serves the same reads, but
its download speed was too low to use. The script skips runs that are already complete.

### 03. Raw read quality control

```bash
THREADS=8 bash scripts/03_qc_raw.sh
```

FastQC on every raw FASTQ and one MultiQC report per dataset: `results/qc/raw/multiqc_{A,B}.html`.

### 04. Read trimming

```bash
THREADS=8 bash scripts/04_trim.sh
```

Cutadapt with dataset-specific options chosen from the raw FastQC reports, then FastQC and MultiQC
again (`results/qc/trimmed/multiqc_{A,B}.html`, including the Cutadapt logs). Trimmed reads go to
`data/trimmed/{A,B}/`.

| Dataset | Options | Reason |
|---|---|---|
| A | `-u 3 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20` | the first 3 nt of read 1 are the template-switching GGG |
| B | `-u 1 -U 1 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT --nextseq-trim=20` | the first base of both reads is a non-biological T; Nextera-type adapter; two-colour (NovaSeq) quality trimming |

Both datasets use `--minimum-length 20`.

### 05. Mapping and strandedness

```bash
THREADS=8 bash scripts/05_map.sh A
```

Builds the HISAT2 index of `sk36_obelisk.fa` once (`data/ref/hisat2/`) and maps the trimmed reads
of the dataset with `--no-spliced-alignment` (bacteria have no introns). Sorted and indexed BAM
files, unmapped reads included, go to `data/bam/<dataset>/`. For each run,
`results/mapping/<dataset>/` holds the HISAT2 log, `samtools flagstat` and `idxstats`, and the
RSeQC `infer_experiment.py` report. The report is computed on host genes only and skipped for
the media controls, which have almost no reads.

**Strandedness.** In dataset A, 0.965–0.994 of the reads follow `1++,1--,2+-,2-+`: read 1 has
the strand of the RNA (featureCounts `-s 1`, kallisto `--fr-stranded`). The value is set by hand
for each dataset in `06_count.sh` and `07_quant.sh`; dataset B is not set yet.

### 06. Read counts and obelisk metrics

```bash
THREADS=8 bash scripts/06_count.sh A
```

- **featureCounts** (`-p --countReadPairs -s 1 -F SAF`): read pairs per gene in
  `results/counts/<dataset>/featurecounts.tsv`. With the default options only uniquely mapped
  pairs are counted, so pairs from the four identical rRNA operons stay unassigned.
- **`scripts/obelisk_metrics.py`**: `results/obelisk/<dataset>/obelisk_metrics.tsv`, one row
  per run. Each read pair is classified once by its primary alignment (read 1, or read 2 when
  only read 2 is mapped) as unmapped, rRNA, other host, obelisk sense or obelisk antisense. The
  strand is that of the RNA molecule, given by `-s`. rRNA pairs are counted here because
  featureCounts leaves them unassigned.

  | Column | Content |
  |---|---|
  | `pct_mapped` | mapped pairs, % of all pairs |
  | `pct_rrna_mapped` | rRNA pairs, % of mapped pairs |
  | `pct_obelisk_total`, `pct_obelisk_mapped`, `pct_obelisk_non_rrna` | obelisk pairs (both strands), % of all, mapped and non-rRNA mapped pairs |
  | `pct_antisense` | antisense obelisk pairs, % of obelisk pairs |
  | `obelisk_tpm`, `obelisk_rank` | obelisk TPM (sense + antisense, from featureCounts) and its rank among the CDS (1 = most abundant) |

- **`scripts/plot_obelisk_coverage.py`**: `results/obelisk/<dataset>/obelisk_coverage.png`
  (with the data in `.tsv`). Read depth along the obelisk for each RNA strand, per million read
  pairs, with the junction extension folded back onto the 1137 nt circle. There is one panel per
  preparation (ObP1 runs only), with the replicate mean and range on a log scale as in Fig. 5B
  of Zheludev et al. (2024), and the Oblin-1 ORF shaded.
- **MultiQC**: `results/mapping/multiqc_<dataset>.html` (HISAT2, samtools, RSeQC and
  featureCounts).

### 07. kallisto quantification

```bash
THREADS=8 bash scripts/07_quant.sh A
```

Builds the kallisto index of `transcripts.fa` once (`data/ref/kallisto/`) and quantifies the ObP1
and ObN1 runs of the dataset with `--fr-stranded` and no bootstraps; media controls are skipped.
For each run, `results/quant/<dataset>/<run>/` holds `abundance.tsv`, `run_info.json` and the
kallisto log, and `results/quant/multiqc_quant_<dataset>.html` summarises them. The output is plain
text (`--plaintext`), so tximport does not need rhdf5. The obelisk sense and antisense targets are
reverse complements of each other and share all their k-mers: only stranded quantification can
tell them apart.

### 08. Differential expression

```bash
bash scripts/08_de.sh A
```

Runs `scripts/de_deseq2.R` once per library preparation (total RNA, rRNA-depleted and
RNaseR-treated). DESeq2 compares ObP1 with ObN1 (design `~ strain`; a positive log2 fold change
means higher in ObP1; DE when padj < 0.05) on two quantification routes:

- featureCounts read-pair counts (step 06), with `DESeqDataSetFromMatrix`;
- kallisto estimated counts (step 07), imported with tximport and `DESeqDataSetFromTximport`.

In both routes the obelisk sense and antisense counts are added up into one gene,
`Obelisk_000003`, and the 12 rRNA genes are left out, since their level reflects the rRNA
depletion (2337 genes). Outputs in `results/de/<dataset>/<prep>/`:

| File | Content |
|---|---|
| `deseq2_featurecounts.tsv`, `deseq2_kallisto.tsv` | DESeq2 results per gene, sorted by padj |
| `size_factors.tsv` | median-of-ratios and library-size factors per run, with and without the obelisk |
| `summary.tsv` | per route: genes tested, host DE genes, obelisk log2FC and padj, Spearman correlation of the log2FC between routes |
| `pca.png` | PCA of the variance-stabilised counts, host genes only |
| `ma.png`, `volcano.png` | MA and volcano plots; values beyond the range of the host genes are drawn as triangles at the edge |
| `lfc_concordance.png` | log2FC from featureCounts against kallisto, host genes with mean count ≥ 10 in both |

## Status

Work in progress.
