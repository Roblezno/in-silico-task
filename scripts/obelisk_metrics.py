#!/usr/bin/env python3
"""Compute Obelisk-S.s abundance and strand metrics for one dataset.

Read pair classes: every read pair is counted once, from the primary alignment
of read 1 (or of read 2 when only read 2 is mapped), as unmapped, host rRNA,
obelisk sense, obelisk antisense or other host. Pairs that map to several
places (the four identical rRNA operons, the 30 nt repeated at both obelisk
ends) are counted once, at their primary alignment.

Strand: it refers to the RNA molecule. With --strand 1 read 1 has the strand of
the RNA and read 2 the opposite one; with --strand 2 it is the other way round
(same convention as featureCounts -s). Sense is the + strand of the obelisk,
the strand of its Oblin-1 ORF.

Expression rank: from the featureCounts matrix, genes are ranked by TPM among
the protein-coding (CDS) genes plus the obelisk (sense + antisense counts).

Output: one row per run with the sample sheet columns and the metrics.
"""

import argparse
from pathlib import Path

import pandas as pd
import pysam

OBELISK_BIOTYPE = "obelisk"
RRNA_BIOTYPE = "rRNA"
CODING_BIOTYPE = "CDS"
PAIR_CLASSES = ["unmapped", "rrna", "other_host", "obelisk_sense", "obelisk_antisense"]


def rrna_intervals(genes: pd.DataFrame) -> dict[str, list[tuple[int, int]]]:
    """Return the rRNA gene intervals per chromosome (0-based, half-open)."""
    intervals: dict[str, list[tuple[int, int]]] = {}
    for gene in genes[genes["biotype"] == RRNA_BIOTYPE].itertuples():
        intervals.setdefault(gene.chrom, []).append((gene.start - 1, gene.end))
    return intervals


def overlaps_any(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    """Return True if [start, end) overlaps any of the intervals."""
    return any(start < i_end and i_start < end for i_start, i_end in intervals)


def rna_is_forward(read: pysam.AlignedSegment, strand: int) -> bool:
    """Return True if the RNA molecule the read comes from is on the + strand."""
    read_has_rna_strand = read.is_read1 == (strand == 1)
    return read.is_reverse != read_has_rna_strand


def classify_read(
    read: pysam.AlignedSegment,
    obelisk_chrom: str,
    rrna: dict[str, list[tuple[int, int]]],
    strand: int,
) -> str:
    """Return the class of the pair represented by this mapped read."""
    if read.reference_name == obelisk_chrom:
        return "obelisk_sense" if rna_is_forward(read, strand) else "obelisk_antisense"
    if overlaps_any(read.reference_start, read.reference_end, rrna.get(read.reference_name, [])):
        return "rrna"
    return "other_host"


def count_pair_classes(
    bam_path: Path,
    obelisk_chrom: str,
    rrna: dict[str, list[tuple[int, int]]],
    strand: int,
    threads: int,
) -> dict[str, int]:
    """Count the read pairs of one BAM file in each class (each pair once)."""
    counts = dict.fromkeys(PAIR_CLASSES, 0)
    with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_read1:
                if read.is_unmapped:
                    # Both unmapped: count here. Only read 2 mapped: count it from read 2.
                    if read.mate_is_unmapped:
                        counts["unmapped"] += 1
                    continue
            elif read.is_unmapped or not read.mate_is_unmapped:
                # Read 2 represents the pair only when read 1 is unmapped.
                continue
            counts[classify_read(read, obelisk_chrom, rrna, strand)] += 1
    return counts


def pair_metrics(pairs: pd.DataFrame) -> pd.DataFrame:
    """Add totals and percentages to the per-run pair class counts."""
    total = pairs[PAIR_CLASSES].sum(axis=1)
    mapped = total - pairs["unmapped"]
    obelisk = pairs["obelisk_sense"] + pairs["obelisk_antisense"]
    metrics = pd.DataFrame({
        "run": pairs["run"],
        "total_pairs": total,
        "mapped_pairs": mapped,
        "pct_mapped": 100 * mapped / total,
        "rrna_pairs": pairs["rrna"],
        "other_host_pairs": pairs["other_host"],
        "obelisk_sense_pairs": pairs["obelisk_sense"],
        "obelisk_antisense_pairs": pairs["obelisk_antisense"],
        "obelisk_pairs": obelisk,
        "pct_rrna_mapped": 100 * pairs["rrna"] / mapped,
        "pct_obelisk_total": 100 * obelisk / total,
        "pct_obelisk_mapped": 100 * obelisk / mapped,
        "pct_obelisk_non_rrna": 100 * obelisk / (mapped - pairs["rrna"]),
        "pct_antisense": 100 * pairs["obelisk_antisense"] / obelisk,
    })
    return metrics


def read_featurecounts(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Return the featureCounts matrix (genes x runs) and the gene lengths."""
    table = pd.read_csv(path, sep="\t", comment="#", index_col="Geneid")
    lengths = table["Length"]
    counts = table.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    counts.columns = [Path(column).stem for column in counts.columns]
    return counts, lengths


def obelisk_expression(counts: pd.DataFrame, lengths: pd.Series, genes: pd.DataFrame) -> pd.DataFrame:
    """Return the obelisk TPM (sense + antisense) and its rank among CDS genes.

    The TPM is missing for runs without assigned reads, and the rank is missing
    when the obelisk is not detected (a zero would only be tied with every
    other unexpressed gene).
    """
    rate = counts.div(lengths, axis=0)
    tpm = rate / rate.sum() * 1e6
    biotype = genes.set_index("gene_id")["biotype"]
    obelisk_ids = biotype.index[biotype == OBELISK_BIOTYPE]
    cds_ids = biotype.index[biotype == CODING_BIOTYPE]

    obelisk_tpm = tpm.loc[obelisk_ids].sum(min_count=1)
    candidates = pd.concat([tpm.loc[cds_ids], obelisk_tpm.to_frame("obelisk").T])
    rank = candidates.rank(ascending=False, method="min").loc["obelisk"].where(obelisk_tpm > 0)
    return pd.DataFrame({
        "obelisk_tpm": obelisk_tpm,
        "obelisk_rank": rank,
        "ranked_genes": len(candidates),
    })


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", type=Path, required=True, help="sample sheet (config/samples_<dataset>.tsv)")
    parser.add_argument("--bamdir", type=Path, required=True, help="directory with <run>.bam files")
    parser.add_argument("--genes", type=Path, required=True, help="gene table (data/ref/genes.tsv)")
    parser.add_argument("--counts", type=Path, required=True, help="featureCounts matrix")
    parser.add_argument("--strand", type=int, choices=[1, 2], required=True,
                        help="read that has the RNA strand (featureCounts -s convention)")
    parser.add_argument("--threads", type=int, default=1, help="BAM decompression threads (default: 1)")
    parser.add_argument("--output", type=Path, required=True, help="output table (TSV)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t").drop(columns="read_count")
    genes = pd.read_csv(args.genes, sep="\t")
    obelisk_chrom = genes.loc[genes["biotype"] == OBELISK_BIOTYPE, "chrom"].iloc[0]
    rrna = rrna_intervals(genes)

    rows = []
    for run in samples["run"]:
        print(f"Classifying read pairs of {run}")
        counts = count_pair_classes(args.bamdir / f"{run}.bam", obelisk_chrom, rrna, args.strand, args.threads)
        rows.append({"run": run, **counts})
    metrics = pair_metrics(pd.DataFrame(rows))

    counts, lengths = read_featurecounts(args.counts)
    expression = obelisk_expression(counts, lengths, genes)

    table = samples.merge(metrics, on="run").merge(expression, left_on="run", right_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.round(3).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
