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
config/    sample sheets and parameters
scripts/   analysis steps, numbered in execution order
report/    Quarto report
data/      downloaded data and references (not versioned)
results/   intermediate and final results (not versioned)
```

## Status

Work in progress.
