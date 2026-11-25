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
THREADS=8 bash scripts/05_map.sh A    # or B
```

Builds the HISAT2 index of `sk36_obelisk.fa` once (`data/ref/hisat2/`) and maps the trimmed reads
of the dataset with `--no-spliced-alignment` (bacteria have no introns). Sorted and indexed BAM
files, unmapped reads included, go to `data/bam/<dataset>/`. For each run,
`results/mapping/<dataset>/` holds the HISAT2 log, `samtools flagstat` and `idxstats`, and the
RSeQC `infer_experiment.py` report. The report is computed on host genes only and skipped for
the media controls, which have almost no reads.

**Strandedness.** The value is set by hand for each dataset in `06_count.sh` and `07_quant.sh`
after reading the `infer_experiment.py` reports:

| Dataset | Reads following the pattern | Read 1 | featureCounts | kallisto |
|---|---|---|---|---|
| A | 0.965–0.994 `1++,1--,2+-,2-+` | sense to the RNA | `-s 1` | `--fr-stranded` |
| B | 0.964–0.997 `1+-,1-+,2++,2--` | antisense to the RNA (dUTP library) | `-s 2` | `--rf-stranded` |

### 06. Read counts and obelisk metrics

```bash
THREADS=8 bash scripts/06_count.sh A    # or B
```

- **featureCounts** (`-p --countReadPairs -F SAF`, with `-s 1` in A and `-s 2` in B): read pairs per gene in
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
  preparation in A (ObP1 runs only) and one per oxygen condition in B (all runs), with the
  replicate mean and range on a log scale as in Fig. 5B of Zheludev et al. (2024), and the
  Oblin-1 ORF shaded.
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

Dataset A only: dataset B is analysed with featureCounts in step 08 and quantified with kallisto
only in the mapper benchmark (step 10, one subsampled run), so `07_quant.sh B` stops with an error.

### 08. Differential expression

```bash
bash scripts/08_de.sh A    # or B
```

In both datasets the obelisk sense and antisense counts are added up into one gene,
`Obelisk_000003`, and the 12 rRNA genes are left out, since their level reflects the rRNA
depletion (2337 genes). A gene is DE when padj < 0.05.

**Dataset A: ObP1 against ObN1.** Runs `scripts/de_deseq2.R` once per library preparation (total
RNA, rRNA-depleted and RNaseR-treated). DESeq2 compares ObP1 with ObN1 (design `~ strain`; a
positive log2 fold change means higher in ObP1) on two quantification routes:

- featureCounts read-pair counts (step 06), with `DESeqDataSetFromMatrix`;
- kallisto estimated counts (step 07), imported with tximport and `DESeqDataSetFromTximport`.

Outputs in `results/de/A/<prep>/`:

| File | Content |
|---|---|
| `deseq2_featurecounts.tsv`, `deseq2_kallisto.tsv` | DESeq2 results per gene, sorted by padj |
| `size_factors.tsv` | median-of-ratios and library-size factors per run, with and without the obelisk |
| `summary.tsv` | per route: genes tested, host DE genes, obelisk log2FC and padj, Spearman correlation of the log2FC between routes |
| `pca.png` | PCA of the variance-stabilised counts, host genes only |
| `ma.png`, `volcano.png` | MA and volcano plots; values beyond the range of the host genes are drawn as triangles at the edge |
| `lfc_concordance.png` | log2FC from featureCounts against kallisto, host genes with mean count ≥ 10 in both |

**Dataset B: culture conditions.** Runs `scripts/de_conditions.R` on the featureCounts counts with
design `~ oxygen + sugar` (no interaction, too few replicates). Two contrasts come from the same
fit: oxic against anoxic and sucrose against glucose (a positive log2 fold change means higher in
oxic or in sucrose). It runs twice:

- `results/de/B/all/`: the 10 runs;
- `results/de/B/no_oxic_rep3/`: without `glucose_oxic_3` and `sucrose_oxic_3` (`--exclude`), a
  sensitivity run because these two runs stand out (22–24 % rRNA and 14–16 % obelisk pairs) and
  may carry a batch effect.

| File | Content |
|---|---|
| `deseq2_oxygen.tsv`, `deseq2_sugar.tsv` | DESeq2 results per gene and contrast, sorted by padj |
| `summary.tsv` | per contrast: samples, genes tested, host DE genes (up and down), obelisk log2FC, SE and padj |
| `obelisk_counts.tsv`, `obelisk_counts.png` | raw and normalised obelisk counts per run, by condition |
| `pca.png` | PCA of the variance-stabilised counts, host genes only; colour = oxygen, shape = sugar |
| `ma.png`, `volcano.png` | MA and volcano plots, one panel per contrast; the volcano has one y axis per panel |

### 09. Comparison with the published values

```bash
bash scripts/09_compare.sh
```

Runs `scripts/compare_studies.py`, which puts the obelisk metrics of step 06 next to the values
published for each dataset. The published values were copied by hand into `resources/`, with
their source:

- `zheludev2024_values.tsv`: Zheludev et al. (2024), mean ± SD per preparation of the obelisk and
  rRNA percentages (Fig. 4D) and of the antisense percentage (Fig. 5B);
- `maddamsetti2025_values.tsv`: Maddamsetti & You (2025), per run, the obelisk rank (Table 1) and
  the mapped and obelisk percentages (Table 2).

Outputs in `results/comparison/`:

| File | Content |
|---|---|
| `obelisk_summary.tsv` | one row per group (A: ObP1 per preparation; B: per oxygen condition), our values as mean ± SD next to the published ones, formatted for the report |
| `b_vs_maddamsetti.tsv` | one row per run of dataset B: our and the published percentages, their ratio and the obelisk ranks |

The published percentages are not defined exactly like ours (reads or k-mers instead of read
pairs, other mappers and filters): the `note` column of the summary says how each one was obtained.

### 10. Mapper benchmark

```bash
THREADS=8 bash scripts/10_benchmark.sh A    # or B
```

Compares how four obelisk reference variants and five mapping configurations handle the circular
obelisk, on one sample per dataset: A uses SRR29668056 (ObP1-A_totalRNA) whole (0.70 M pairs), B
uses SRR23591557 (sucrose, oxic, replicate 1) subsampled to 2 M pairs with `seqtk sample -s 42`
(the same seed on both mates keeps the pairs). There are no simulated reads, so the benchmark
measures agreement between tools, not accuracy.

| Variant | Obelisk contig | Note |
|---|---|---|
| R0 | 1137 nt | linear: no read can cross the circle junction |
| R1 | + 30 nt | the reference of steps 01–09 |
| R2 | + read length − 1 (A: 77, B: 98) | every read crossing the junction fits in one piece |
| R3 | + 1137 nt | dimer: every obelisk read maps to two places |

The variants are built with `build_references.py --extension` in `data/ref/benchmark/<dataset>/`,
each with its BWA, Bowtie2, HISAT2 and kallisto index. The tools run with default options except:
`bwa mem`, `bowtie2` (end-to-end), `bowtie2 --local`, `hisat2 --no-spliced-alignment` and
`kallisto quant` with the strand of the dataset (`--fr-stranded` in A, `--rf-stranded` in B).
GNU time (`time -v`) records the run time and peak memory of the mapper process only. References,
indexes and mappings are made once; a finished mapping leaves a `.done` marker in
`results/benchmark/<dataset>/logs/` (delete it to re-map). Counting, metrics and figures are
always re-run.

- **featureCounts** per variant on the four BAM files (`-M --primary`, strand as in step 06):
  every primary alignment is counted, so pairs that map to repeated obelisk sequence (R2, R3)
  count once. Supplementary alignments are counted too: bwa splits reads that cross the junction
  in R0, which raises its obelisk TPM.
- **`scripts/obelisk_metrics.py`**, unchanged, with one sheet row per aligner (`run` =
  `<ref>.<tool>`): the obelisk metrics of step 06.
- **`scripts/benchmark_metrics.py`**: `results/benchmark/<dataset>/benchmark_summary.tsv`, one row
  per variant × tool, and `host_spearman_R1.tsv`.

  | Column | Content |
  |---|---|
  | `pct_obelisk_total`, `pct_antisense`, `obelisk_tpm`, `obelisk_rank` | as in step 06; for kallisto from `est_counts`, and "mapped" means pseudoaligned |
  | `junction_reads`, `pct_junction_reads` | primary obelisk reads with aligned bases on both sides of the circle junction |
  | `end_clipped_reads` | primary obelisk reads soft-clipped at a contig end (in R0, reads cut by the junction) |
  | `supplementary_alignments` | supplementary alignments on the obelisk contig |
  | `pct_mapq0`, `median_mapq` | MAPQ of the primary obelisk reads (scales differ between aligners) |
  | `pct_softclipped`, `mean_softclip_nt` | primary obelisk reads with soft clipping, and soft-clipped nt per read |
  | `wall_s`, `user_s`, `sys_s`, `cpu_pct`, `max_rss_mb` | GNU time report of the mapper |

  The alignment columns are empty for kallisto. `host_spearman_R1.tsv` is the Spearman
  correlation of the host CDS counts between the five tools with variant R1.
- **`scripts/plot_benchmark.py`**: four figures in `results/benchmark/<dataset>/figures/`:
  `obelisk_abundance.png` (obelisk and antisense shares), `obelisk_junction.png` (junction and
  end-clipped reads), `obelisk_mapq_softclip.png` (median MAPQ and soft clipping) and
  `host_spearman_resources.png` (host correlation, wall time and peak memory with R1).

## Status

Work in progress.
