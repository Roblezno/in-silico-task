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
| `annotation.saf` | host genes + obelisk | featureCounts |
| `host_genes.bed` | BED12, host genes only | RSeQC `infer_experiment.py` |
| `transcripts.fa` | one sequence per gene on its own strand + obelisk | kallisto |
| `genes.tsv` | gene ID, biotype, name, product, coordinates, length | metrics, report |

**Obelisk-S.s sequence.** `resources/obelisk_000003.fasta` holds Obelisk_000003 (1137 nt),
extracted once from the supplementary table `mmc2.xlsx` of Zheludev et al. (2024) with
`scripts/extract_obelisk.py` (row `Obelisk_000003_000001_000001`). Its Oblin-1 ORF lies on
the + strand (nt 51–743), which defines the obelisk *sense* strand. The obelisk is circular:
its first 30 nt are appended to its end so that reads crossing the circle junction can align,
as in the original article.

## Status

Work in progress.
