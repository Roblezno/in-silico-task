#!/usr/bin/env python3
"""Summarize the mapper benchmark of one dataset (4 reference variants x 5 tools).

Variants R0-R3 are the host genome plus the obelisk extended by 0, 30, a read
length and 1137 nt. Tools are four aligners (bwa mem, bowtie2 end-to-end,
bowtie2 --local, HISAT2) and kallisto.

Obelisk abundance: the same metrics as the main analysis. For the aligners
they are read from the table of scripts/obelisk_metrics.py (pairs counted
once, at their primary alignment). For kallisto the sense and antisense
obelisk counts are the est_counts of the two obelisk targets (fractional
counts), and "mapped" means pseudoaligned: pct_mapped is n_pseudoaligned
over n_processed and pct_obelisk_mapped is the obelisk share of the
pseudoaligned pairs. The obelisk TPM and rank are computed with
obelisk_expression from est_counts and eff_length, which reproduces the TPM
of kallisto.

Obelisk alignment metrics (aligners only, empty for kallisto): from the BAM,
on the obelisk contig only. Counts are per read (each mate once), not per
pair.
- obelisk_alignments: primary mapped reads (not secondary, not supplementary).
- junction_reads: primary reads with aligned bases on both sides of the circle
  junction, i.e. an aligned block covers the 0-based positions L-1 and L, with
  L the length of the unextended obelisk. Always 0 in R0, whose contig ends at L.
- end_clipped_reads: primary reads soft-clipped at a contig end (leading clip
  starting at position 0, or trailing clip ending at the contig end). In R0
  these are reads cut by the junction.
- supplementary_alignments: supplementary alignments (flag 0x800) on the
  contig (bwa splits reads that cross the junction in R0).
- pct_mapq0, median_mapq: MAPQ 0 share and median MAPQ of the primary reads.
- pct_softclipped, mean_softclip_nt: share of primary reads with any soft clip,
  and mean soft-clipped nt per primary read (clipped or not).

Resources (all tools): wall time, user and system CPU time, CPU use and peak
memory from the GNU time -v report of the mapper process.

Host correlation: Spearman correlation matrix, between the five tools, of the
host gene counts of variant R1 over the protein-coding (CDS) genes. Ties (many
genes with zero counts) are handled by average ranks.

Outputs (in --resultdir):
  benchmark_summary.tsv   one row per variant x tool
  host_spearman_R1.tsv    5 x 5 Spearman correlation matrix
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import pysam

from obelisk_metrics import CODING_BIOTYPE, OBELISK_BIOTYPE, obelisk_expression, read_featurecounts

VARIANTS = ["R0", "R1", "R2", "R3"]
ALIGNERS = ["bwa_mem", "bowtie2_e2e", "bowtie2_local", "hisat2"]
KALLISTO = "kallisto"
SPEARMAN_VARIANT = "R1"
SOFT_CLIP = 4  # CIGAR operation S in pysam.AlignedSegment.cigartuples

ABUNDANCE_COLUMNS = [
    "pct_mapped", "obelisk_sense_pairs", "obelisk_antisense_pairs", "obelisk_pairs",
    "pct_obelisk_total", "pct_obelisk_mapped", "pct_antisense",
    "obelisk_tpm", "obelisk_rank", "ranked_genes",
]
ALIGNMENT_COLUMNS = [
    "obelisk_alignments", "junction_reads", "pct_junction_reads", "end_clipped_reads",
    "supplementary_alignments", "pct_mapq0", "median_mapq", "pct_softclipped", "mean_softclip_nt",
]
RESOURCE_COLUMNS = ["wall_s", "user_s", "sys_s", "cpu_pct", "max_rss_mb"]
SUMMARY_COLUMNS = (
    ["ref", "tool", "extension", "input_pairs"] + ABUNDANCE_COLUMNS + ALIGNMENT_COLUMNS + RESOURCE_COLUMNS
)
# Columns that hold whole numbers (nullable, so kallisto rows stay empty).
INTEGER_COLUMNS = [
    "obelisk_rank", "ranked_genes", "obelisk_alignments", "junction_reads",
    "end_clipped_reads", "supplementary_alignments",
]

# Lines of the GNU time -v report that are used.
TIME_FIELDS = {
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall_s",
    "User time (seconds)": "user_s",
    "System time (seconds)": "sys_s",
    "Percent of CPU this job got": "cpu_pct",
    "Maximum resident set size (kbytes)": "max_rss_mb",
}
READ_COLUMNS = ["mapq", "softclip_nt", "junction", "end_clipped"]


def percent(part: float, whole: float) -> float:
    """Return 100 * part / whole, or NaN when whole is zero."""
    return 100 * part / whole if whole else float("nan")


def fasta_length(path: Path) -> int:
    """Return the number of nucleotides of a single-record FASTA file."""
    with open(path) as handle:
        return sum(len(line.strip()) for line in handle if not line.startswith(">"))


def obelisk_contig(genes: pd.DataFrame) -> tuple[str, int]:
    """Return the name and length of the obelisk contig of a variant."""
    row = genes[genes["biotype"] == OBELISK_BIOTYPE].iloc[0]
    return row["chrom"], int(row["length"])


# --- Time report ---------------------------------------------------------


def parse_elapsed(text: str) -> float:
    """Convert a GNU time elapsed field (h:mm:ss or m:ss.ss) to seconds."""
    seconds = 0.0
    for part in text.split(":"):
        seconds = seconds * 60 + float(part)
    return seconds


def parse_time_report(text: str) -> dict[str, float]:
    """Return wall, user and system time (s), CPU use (%) and peak memory (MB)."""
    raw = {}
    for line in text.splitlines():
        key, _, value = line.strip().partition(": ")
        if key in TIME_FIELDS:
            raw[TIME_FIELDS[key]] = value
    missing = set(RESOURCE_COLUMNS) - raw.keys()
    if missing:
        raise ValueError(f"fields missing in the time report: {sorted(missing)}")
    return {
        "wall_s": parse_elapsed(raw["wall_s"]),
        "user_s": float(raw["user_s"]),
        "sys_s": float(raw["sys_s"]),
        "cpu_pct": float(raw["cpu_pct"].rstrip("%")),
        "max_rss_mb": float(raw["max_rss_mb"]) / 1024,
    }


# --- Obelisk alignments of one BAM ---------------------------------------


def spans_junction(blocks: list[tuple[int, int]], junction: int) -> bool:
    """Return True if aligned blocks cover the 0-based positions junction-1 and junction."""
    left = any(start <= junction - 1 < end for start, end in blocks)
    right = any(start <= junction < end for start, end in blocks)
    return left and right


def read_features(read: pysam.AlignedSegment, contig_length: int, junction: int) -> dict:
    """Return the features of one primary mapped read on the obelisk contig."""
    cigar = read.cigartuples
    leading = cigar[0][1] if cigar[0][0] == SOFT_CLIP else 0
    trailing = cigar[-1][1] if cigar[-1][0] == SOFT_CLIP else 0
    return {
        "mapq": read.mapping_quality,
        "softclip_nt": leading + trailing,
        "junction": spans_junction(read.get_blocks(), junction),
        "end_clipped": (leading > 0 and read.reference_start == 0)
        or (trailing > 0 and read.reference_end == contig_length),
    }


def read_obelisk_alignments(
    bam_path: Path, contig: str, junction: int, threads: int
) -> tuple[pd.DataFrame, int]:
    """Return the features of the primary reads on the contig and the supplementary count."""
    records = []
    supplementary = 0
    with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
        contig_length = bam.get_reference_length(contig)
        for read in bam.fetch(contig):
            # An unmapped mate is placed at the position of its mapped mate.
            if read.is_unmapped:
                continue
            if read.is_supplementary:
                supplementary += 1
            elif not read.is_secondary:
                records.append(read_features(read, contig_length, junction))
    reads = pd.DataFrame(records, columns=READ_COLUMNS)
    reads = reads.astype({"mapq": int, "softclip_nt": int, "junction": bool, "end_clipped": bool})
    return reads, supplementary


def summarize_reads(reads: pd.DataFrame, supplementary: int) -> dict[str, float]:
    """Return the obelisk alignment metrics from the features of the primary reads."""
    return {
        "obelisk_alignments": len(reads),
        "junction_reads": int(reads["junction"].sum()),
        "pct_junction_reads": 100 * reads["junction"].mean(),
        "end_clipped_reads": int(reads["end_clipped"].sum()),
        "supplementary_alignments": supplementary,
        "pct_mapq0": 100 * (reads["mapq"] == 0).mean(),
        "median_mapq": reads["mapq"].median(),
        "pct_softclipped": 100 * (reads["softclip_nt"] > 0).mean(),
        "mean_softclip_nt": reads["softclip_nt"].mean(),
    }


# --- kallisto ------------------------------------------------------------


def read_kallisto(directory: Path) -> tuple[pd.DataFrame, int, int]:
    """Return the abundance table (by target) and n_processed, n_pseudoaligned."""
    abundance = pd.read_csv(directory / "abundance.tsv", sep="\t", index_col="target_id")
    info = json.loads((directory / "run_info.json").read_text())
    return abundance, info["n_processed"], info["n_pseudoaligned"]


def kallisto_abundance_metrics(
    abundance: pd.DataFrame, n_processed: int, n_pseudoaligned: int, genes: pd.DataFrame, name: str
) -> dict[str, float]:
    """Return the obelisk abundance metrics of a kallisto run (mapped = pseudoaligned)."""
    obelisk_ids = genes[genes["biotype"] == OBELISK_BIOTYPE].set_index("strand")["gene_id"]
    sense = abundance.loc[obelisk_ids["+"], "est_counts"]
    antisense = abundance.loc[obelisk_ids["-"], "est_counts"]
    obelisk = sense + antisense
    counts = abundance[["est_counts"]].rename(columns={"est_counts": name})
    expression = obelisk_expression(counts, abundance["eff_length"], genes).loc[name]
    return {
        "pct_mapped": percent(n_pseudoaligned, n_processed),
        "obelisk_sense_pairs": sense,
        "obelisk_antisense_pairs": antisense,
        "obelisk_pairs": obelisk,
        "pct_obelisk_total": percent(obelisk, n_processed),
        "pct_obelisk_mapped": percent(obelisk, n_pseudoaligned),
        "pct_antisense": percent(antisense, obelisk),
        **expression.to_dict(),
    }


# --- Summary table -------------------------------------------------------


def aligner_row(
    run: str, obelisk_metrics: pd.DataFrame, bam_path: Path, contig: str, junction: int, threads: int
) -> dict[str, float]:
    """Return the abundance and obelisk alignment metrics of one aligner."""
    print(f"Reading obelisk alignments of {bam_path.name}")
    reads, supplementary = read_obelisk_alignments(bam_path, contig, junction, threads)
    return {
        **obelisk_metrics.loc[run, ABUNDANCE_COLUMNS].to_dict(),
        **summarize_reads(reads, supplementary),
    }


def kallisto_row(ref: str, kallisto_dir: Path, genes: pd.DataFrame) -> dict[str, float]:
    """Return the abundance metrics of kallisto; the alignment metrics are missing."""
    abundance, n_processed, n_pseudoaligned = read_kallisto(kallisto_dir)
    return {
        **kallisto_abundance_metrics(abundance, n_processed, n_pseudoaligned, genes, f"{ref}.{KALLISTO}"),
        **dict.fromkeys(ALIGNMENT_COLUMNS, float("nan")),
    }


def summarize_variant(args: argparse.Namespace, ref: str, unextended: int, input_pairs: int) -> list[dict]:
    """Return the summary rows (one per tool) of one reference variant."""
    genes = pd.read_csv(args.refdir / ref / "genes.tsv", sep="\t")
    contig, contig_length = obelisk_contig(genes)
    metrics_path = args.resultdir / "metrics" / f"{ref}.obelisk_metrics.tsv"
    obelisk_metrics = pd.read_csv(metrics_path, sep="\t", index_col="run")
    rows = []
    for tool in ALIGNERS + [KALLISTO]:
        row = {"ref": ref, "tool": tool, "extension": contig_length - unextended, "input_pairs": input_pairs}
        if tool == KALLISTO:
            row.update(kallisto_row(ref, args.resultdir / KALLISTO / ref, genes))
        else:
            run = f"{ref}.{tool}"
            bam_path = args.bamdir / f"{run}.bam"
            row.update(aligner_row(run, obelisk_metrics, bam_path, contig, unextended, args.threads))
        report = args.resultdir / "logs" / f"{ref}.{tool}.time.txt"
        row.update(parse_time_report(report.read_text()))
        rows.append(row)
    return rows


# --- Host gene correlation -----------------------------------------------


def host_spearman(
    counts: pd.DataFrame, kallisto_counts: pd.Series, genes: pd.DataFrame, ref: str
) -> pd.DataFrame:
    """Return the Spearman correlation matrix of the host CDS counts between tools.

    counts has the featureCounts columns named <ref>.<tool>; kallisto_counts is
    the kallisto est_counts by target.
    """
    table = counts.rename(columns=lambda column: column.removeprefix(f"{ref}."))
    table[KALLISTO] = kallisto_counts
    cds_ids = genes.loc[genes["biotype"] == CODING_BIOTYPE, "gene_id"]
    return table.loc[cds_ids, ALIGNERS + [KALLISTO]].corr(method="spearman")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bamdir", type=Path, required=True, help="directory with <ref>.<tool>.bam files")
    parser.add_argument("--refdir", type=Path, required=True, help="data/ref/benchmark/<dataset>")
    parser.add_argument("--resultdir", type=Path, required=True, help="results/benchmark/<dataset>")
    parser.add_argument("--obelisk", type=Path, required=True,
                        help="unextended obelisk FASTA (gives the circle junction position)")
    parser.add_argument("--threads", type=int, default=1, help="BAM decompression threads (default: 1)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    unextended = fasta_length(args.obelisk)
    input_pairs = int(pd.read_csv(args.resultdir / "input_pairs.tsv", sep="\t")["pairs"].iloc[0])

    rows = []
    for ref in VARIANTS:
        rows.extend(summarize_variant(args, ref, unextended, input_pairs))
    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS).round(3)
    summary = summary.astype({column: "Int64" for column in INTEGER_COLUMNS})
    summary.to_csv(args.resultdir / "benchmark_summary.tsv", sep="\t", index=False)

    variant = SPEARMAN_VARIANT
    counts, _ = read_featurecounts(args.resultdir / "counts" / f"{variant}.featurecounts.tsv")
    abundance, _, _ = read_kallisto(args.resultdir / KALLISTO / variant)
    genes = pd.read_csv(args.refdir / variant / "genes.tsv", sep="\t")
    correlation = host_spearman(counts, abundance["est_counts"], genes, variant)
    correlation.round(3).to_csv(args.resultdir / f"host_spearman_{variant}.tsv", sep="\t", index_label="tool")


if __name__ == "__main__":
    main()
