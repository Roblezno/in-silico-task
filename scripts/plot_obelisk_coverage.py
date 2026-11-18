#!/usr/bin/env python3
"""Plot the read coverage along Obelisk-S.s, split by RNA strand.

For each selected run, every read with a primary alignment on the obelisk
contig adds one to each reference position it covers (both reads of a pair
count, as in a genome browser). The strand is that of the RNA molecule, with
the same convention as obelisk_metrics.py. The nt appended to the contig to
cover the circle junction are folded back onto the first positions, so the
coverage is given on the circular genome (1137 nt).

The depth is scaled to depth per million read pairs of the run. The figure has
one panel per group (e.g. library preparation): for each strand the replicate
mean is drawn as a line and the replicate range as a band, on a log scale as in
Fig. 5B of Zheludev et al. (2024), with the Oblin-1 ORF shaded. Zero depth
cannot be drawn on a log scale and is left blank.

Outputs:
    <prefix>.tsv   coverage per run, strand and position (table view of the figure)
    <prefix>.png   figure
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pysam  # noqa: E402
from Bio import SeqIO  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402

from obelisk_metrics import rna_is_forward  # noqa: E402

# Oblin-1 ORF on the + strand, 1-based inclusive (checked in phase F1).
OBLIN1_ORF = (51, 743)
STRAND_COLORS = {"antisense": "#2a78d6", "sense": "#eb6834"}
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
ORF_FILL = "#f0efec"


def strand_depth(bam_path: Path, contig: str, strand: int, threads: int) -> dict[str, np.ndarray]:
    """Return the per-position read depth of each RNA strand along the contig."""
    with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
        length = bam.get_reference_length(contig)
        depth = {"antisense": np.zeros(length, dtype=np.int64), "sense": np.zeros(length, dtype=np.int64)}
        for read in bam.fetch(contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            key = "sense" if rna_is_forward(read, strand) else "antisense"
            for start, end in read.get_blocks():
                depth[key][start:end] += 1
    return depth


def fold_circle(depth: np.ndarray, circle_length: int) -> np.ndarray:
    """Add the depth of the appended junction copy back onto the circle start."""
    folded = depth[:circle_length].copy()
    junction_copy = depth[circle_length:]
    folded[: len(junction_copy)] += junction_copy
    return folded


def select_runs(samples: pd.DataFrame, keep: str | None, group_by: str, groups: list[str]) -> pd.DataFrame:
    """Return the sample sheet rows that pass the --keep filter and belong to the groups."""
    if keep:
        column, value = keep.split("=", 1)
        samples = samples[samples[column] == value]
    return samples[samples[group_by].isin(groups)]


def coverage_table(
    runs: pd.Series, bamdir: Path, contig: str, circle_length: int,
    total_pairs: pd.Series, strand: int, threads: int,
) -> pd.DataFrame:
    """Return the long-format coverage table for the selected runs."""
    positions = np.arange(1, circle_length + 1)
    tables = []
    for run in runs:
        print(f"Obelisk coverage of {run}")
        depth = strand_depth(bamdir / f"{run}.bam", contig, strand, threads)
        for strand_name, values in depth.items():
            folded = fold_circle(values, circle_length)
            tables.append(pd.DataFrame({
                "run": run,
                "strand": strand_name,
                "position": positions,
                "depth": folded,
                "depth_per_million_pairs": folded * 1e6 / total_pairs[run],
            }))
    return pd.concat(tables, ignore_index=True)


def style_axis(ax: Axes) -> None:
    """Recessive axes: no top/right spines, hairline grid, muted ticks."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(which="both", colors=INK_MUTED, labelcolor=INK_SECONDARY, labelsize=8)
    ax.grid(axis="y", which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def plot_coverage(
    coverage: pd.DataFrame, runs: pd.DataFrame, pct_antisense: pd.Series,
    group_by: str, groups: list[str], circle_length: int, path: Path,
) -> None:
    """Draw one coverage panel per group and save the figure."""
    group_of_run = runs.set_index("run")[group_by]
    coverage = coverage.assign(group=coverage["run"].map(group_of_run))
    fig, axes = plt.subplots(len(groups), 1, figsize=(7.5, 2.3 * len(groups)), sharex=True, layout="constrained")
    axes = np.atleast_1d(axes)

    for ax, group in zip(axes, groups):
        data = coverage[coverage["group"] == group]
        group_runs = group_of_run.index[group_of_run == group]
        ax.axvspan(*OBLIN1_ORF, color=ORF_FILL, linewidth=0, zorder=0)
        for strand_name in ("antisense", "sense"):
            wide = data[data["strand"] == strand_name].pivot(
                index="position", columns="run", values="depth_per_million_pairs")
            mean, low, high = (values.where(values > 0) for values in
                               (wide.mean(axis=1), wide.min(axis=1), wide.max(axis=1)))
            color = STRAND_COLORS[strand_name]
            ax.fill_between(wide.index, low, high, color=color, alpha=0.25, linewidth=0)
            ax.plot(wide.index, mean, color=color, linewidth=1.5, label=strand_name)
        ax.set_yscale("log")
        style_axis(ax)
        ax.set_title(
            f"{group}  ·  n = {len(group_runs)}  ·  {pct_antisense[group_runs].mean():.1f}% antisense",
            loc="left", fontsize=9, color=INK_PRIMARY,
        )

    axes[0].set_title(f"shaded: Oblin-1 ORF ({OBLIN1_ORF[0]}–{OBLIN1_ORF[1]})",
                      loc="right", fontsize=8, color=INK_MUTED)
    axes[-1].set_xlim(1, circle_length)
    axes[-1].set_xlabel("Position on obelisk-S.s (nt)", fontsize=9, color=INK_SECONDARY)
    fig.supylabel("Read depth per million read pairs (log scale)", fontsize=9, color=INK_SECONDARY)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper right", ncols=2,
               frameon=False, fontsize=8, labelcolor=INK_SECONDARY)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", type=Path, required=True, help="sample sheet (config/samples_<dataset>.tsv)")
    parser.add_argument("--bamdir", type=Path, required=True, help="directory with <run>.bam files")
    parser.add_argument("--obelisk", type=Path, required=True, help="obelisk FASTA, not extended (circle length)")
    parser.add_argument("--metrics", type=Path, required=True, help="obelisk_metrics.py output (total pairs, % antisense)")
    parser.add_argument("--strand", type=int, choices=[1, 2], required=True,
                        help="read that has the RNA strand (featureCounts -s convention)")
    parser.add_argument("--keep", help="keep only runs with COLUMN=VALUE in the sample sheet (e.g. strain=ObP1)")
    parser.add_argument("--group-by", required=True, help="sample sheet column with one panel per value")
    parser.add_argument("--groups", nargs="+", required=True, help="values of --group-by to plot, in panel order")
    parser.add_argument("--threads", type=int, default=1, help="BAM decompression threads (default: 1)")
    parser.add_argument("--output-prefix", type=Path, required=True, help="prefix for the .tsv and .png outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t")
    runs = select_runs(samples, args.keep, args.group_by, args.groups)
    metrics = pd.read_csv(args.metrics, sep="\t", index_col="run")
    obelisk = SeqIO.read(args.obelisk, "fasta")

    coverage = coverage_table(runs["run"], args.bamdir, obelisk.id, len(obelisk),
                              metrics["total_pairs"], args.strand, args.threads)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(args.output_prefix.with_suffix(".tsv"), sep="\t", index=False, float_format="%.3f")
    plot_coverage(coverage, runs, metrics["pct_antisense"], args.group_by, args.groups,
                  len(obelisk), args.output_prefix.with_suffix(".png"))


if __name__ == "__main__":
    main()
