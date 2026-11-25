#!/usr/bin/env python3
"""Plot the mapper benchmark of one dataset (4 reference variants x 5 tools).

Reads the tables of scripts/benchmark_metrics.py from --resultdir. In every
figure a tool keeps its color, marker and position (BWA-MEM, Bowtie2 end-to-end,
Bowtie2 local, HISAT2, kallisto). The line panels show the four variants R0-R3
on an ordinal axis; tools with identical values are shifted by a small fixed
offset so that all markers stay visible.

Outputs (in <resultdir>/figures/):
  obelisk_abundance.png         obelisk share of all pairs and antisense share, 5 tools
  obelisk_junction.png          reads crossing the circle junction and reads clipped
                                at a contig end, aligners only
  obelisk_mapq_softclip.png     median MAPQ and soft-clipped reads, aligners only
  host_spearman_resources.png   Spearman correlation of host gene counts between
                                tools, wall time and peak memory (variant R1)
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from benchmark_metrics import ALIGNERS, KALLISTO, SPEARMAN_VARIANT, VARIANTS  # noqa: E402

TOOLS = ALIGNERS + [KALLISTO]
NAMES = {
    "bwa_mem": "BWA-MEM", "bowtie2_e2e": "Bowtie2 end-to-end", "bowtie2_local": "Bowtie2 local",
    "hisat2": "HISAT2", "kallisto": "kallisto",
}
COLORS = {
    "bwa_mem": "#257de0", "bowtie2_e2e": "#f47230", "bowtie2_local": "#00b894",
    "hisat2": "#f2a600", "kallisto": "#ee7fac",
}
# Secondary encoding: several colors are below 3:1 contrast on white.
MARKERS = {"bwa_mem": "o", "bowtie2_e2e": "s", "bowtie2_local": "D", "hisat2": "^", "kallisto": "v"}
DODGE = 0.07  # horizontal offset between neighboring tools, in variant positions
SPEARMAN_CMAP = LinearSegmentedColormap.from_list(
    "spearman_blue", ["#cce3fd", "#84b9f3", "#338ced", "#1b64be", "#0e3d7e"])

INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

NO_KALLISTO = "kallisto reports no alignments and is not shown."


def style_axis(ax: Axes, grid_axis: str = "y") -> None:
    """Recessive axes: no top/right spines, hairline grid, muted ticks."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(which="both", colors=INK_MUTED, labelcolor=INK_SECONDARY, labelsize=8)
    ax.grid(axis=grid_axis, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def set_panel_title(ax: Axes, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=9, color=INK_PRIMARY)


def add_caption(fig: Figure, lines: list[str]) -> None:
    """Muted note at the bottom left of the figure (constrained layout reserves the space)."""
    if lines:
        fig.supxlabel("\n".join(lines), x=0.01, horizontalalignment="left", fontsize=8, color=INK_MUTED)


def save(fig: Figure, path: Path) -> None:
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


# --- Line panels (figures 1-3) -------------------------------------------


def widen_ylim(ax: Axes, values: pd.Series, min_span: float) -> None:
    """Make the y range at least min_span x the mean value wide, so tiny differences do not look large.

    Shares cannot exceed 100 %, so a range that would pass 100 is shifted down.
    """
    low, high = values.min(), values.max()
    span = min_span * values.mean()
    if high - low >= span:
        return
    middle = (low + high) / 2
    low, high = middle - span / 2, middle + span / 2
    if high > 100:
        low, high = low - (high - 100), 100
    ax.set_ylim(low, high)


def line_panel(
    ax: Axes, summary: pd.DataFrame, column: str, tools: list[str], title: str, from_zero: bool, min_span: float
) -> None:
    """Draw one line per tool across the four variants, dodged so that equal values stay visible."""
    extension = summary.drop_duplicates("ref").set_index("ref")["extension"]
    for tool in tools:
        values = summary[summary["tool"] == tool].set_index("ref").loc[VARIANTS, column]
        # The offset depends on the tool, not on the panel, so it is the same in every figure.
        offset = DODGE * (TOOLS.index(tool) - (len(TOOLS) - 1) / 2)
        ax.plot(np.arange(len(VARIANTS)) + offset, values.to_numpy(), color=COLORS[tool], marker=MARKERS[tool],
                markersize=7, markeredgecolor="white", markeredgewidth=1, linewidth=1.75, label=NAMES[tool])
    ax.set_xticks(range(len(VARIANTS)), [f"{ref}\n+{int(extension[ref])} nt" for ref in VARIANTS])
    ax.set_xlim(-0.5, len(VARIANTS) - 0.5)
    ax.margins(y=0.1)
    if from_zero:
        ax.set_ylim(bottom=0)
    if min_span:
        widen_ylim(ax, summary.loc[summary["tool"].isin(tools), column], min_span)
    style_axis(ax)
    set_panel_title(ax, title)


def line_figure(
    summary: pd.DataFrame, tools: list[str], panels: list[tuple[str, str, bool, float]],
    header: str, notes: list[str], path: Path,
) -> None:
    """Two line panels (column, title, from_zero, min_span) with the figure title and the shared legend above them.

    The figure title is the title of the legend: a suptitle and an "outside
    upper" figure legend are both anchored at the top of the figure and may
    overlap with constrained layout.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), layout="constrained")
    for ax, (column, title, from_zero, min_span) in zip(axes, panels):
        line_panel(ax, summary, column, tools, title, from_zero, min_span)
    handles, labels = axes[0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, loc="outside upper left", ncols=len(tools), frameon=False, fontsize=8,
                        labelcolor=INK_SECONDARY, title=header, title_fontsize=10, alignment="left")
    legend.get_title().set_color(INK_PRIMARY)
    add_caption(fig, notes)
    save(fig, path)


def plot_abundance(summary: pd.DataFrame, title: str, path: Path) -> None:
    """Obelisk share of all pairs and antisense share of the obelisk pairs, all tools.

    The y ranges are at least 20 % (abundance) and 2 % (antisense) of the mean
    value wide, so that differences in the third decimal do not look large.
    """
    panels = [
        ("pct_obelisk_total", "Obelisk share of all read pairs (%)", False, 0.2),
        ("pct_antisense", "Antisense share of obelisk pairs (%)", False, 0.02),
    ]
    line_figure(summary, TOOLS, panels, f"{title}: obelisk abundance by reference variant", [], path)


def plot_junction(summary: pd.DataFrame, title: str, path: Path) -> None:
    """Reads crossing the circle junction and reads clipped at a contig end, aligners only."""
    panels = [
        ("pct_junction_reads", "Obelisk reads crossing the circle junction (%)", True, 0.0),
        ("pct_end_clipped", "Obelisk reads soft-clipped at a contig end (%)", True, 0.0),
    ]
    line_figure(summary, ALIGNERS, panels, f"{title}: obelisk reads at the circle junction", [NO_KALLISTO], path)


def plot_mapq_softclip(summary: pd.DataFrame, title: str, path: Path) -> None:
    """Median MAPQ and soft-clipped reads of the obelisk, aligners only."""
    panels = [
        ("median_mapq", "Median MAPQ of obelisk reads", True, 0.0),
        ("pct_softclipped", "Obelisk reads with soft clipping (%)", True, 0.0),
    ]
    notes = ["MAPQ scales differ between aligners (BWA and HISAT2 max 60, Bowtie2 max 42–44).", NO_KALLISTO]
    line_figure(summary, ALIGNERS, panels, f"{title}: obelisk MAPQ and soft clipping", notes, path)


# --- Host correlation and resources (figure 4) ---------------------------


def draw_heatmap(fig: Figure, ax: Axes, matrix: pd.DataFrame) -> None:
    """Annotated heatmap of the tool-by-tool Spearman matrix with a small colorbar."""
    values = matrix.loc[TOOLS, TOOLS].to_numpy()
    off_diagonal = values[~np.eye(len(TOOLS), dtype=bool)]
    vmin = np.floor(np.round(off_diagonal.min() * 100, 6)) / 100  # e.g. 0.994 -> 0.99
    norm = Normalize(vmin=vmin, vmax=1.0)
    image = ax.imshow(values, cmap=SPEARMAN_CMAP, norm=norm, aspect="auto")
    for row, col in np.ndindex(values.shape):
        # Colors from the middle of the ramp upwards are dark enough for white text.
        color = "white" if norm(values[row, col]) > 0.6 else INK_PRIMARY
        ax.text(col, row, f"{values[row, col]:.3f}", ha="center", va="center", fontsize=8, color=color)
    names = [NAMES[tool] for tool in TOOLS]
    ax.set_xticks(range(len(TOOLS)), names, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(TOOLS)), names)
    ax.set_xticks(np.arange(-0.5, len(TOOLS)), minor=True)  # white gaps between cells
    ax.set_yticks(np.arange(-0.5, len(TOOLS)), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0, labelcolor=INK_SECONDARY, labelsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    set_panel_title(ax, "Spearman correlation of host gene counts (CDS)")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.8, aspect=12, pad=0.03,
                            ticks=np.linspace(vmin, 1.0, 3), format="%.3f")
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(labelsize=8, labelcolor=INK_SECONDARY, color=INK_MUTED, length=2)


def draw_bars(ax: Axes, values: pd.Series, label: str, unit: str, decimals: int) -> None:
    """Horizontal bars from zero, one per tool in the fixed order from top to bottom, with value labels.

    The axis stays linear even when the values span orders of magnitude: bar
    length must be proportional to the value; the labels give the small values.
    """
    positions = np.arange(len(TOOLS))
    bars = ax.barh(positions, values.to_numpy(), height=0.6, color=[COLORS[tool] for tool in TOOLS])
    ax.bar_label(bars, labels=[f"{value:.{decimals}f}" for value in values], padding=3,
                 fontsize=8, color=INK_SECONDARY)
    ax.set_xlim(0, values.max() * 1.2)  # room on the right for the labels
    ax.set_yticks(positions, [NAMES[tool] for tool in TOOLS])
    ax.invert_yaxis()
    style_axis(ax, grid_axis="x")
    ax.tick_params(axis="y", length=0)
    set_panel_title(ax, f"{label} ({unit})")


def plot_host_resources(summary: pd.DataFrame, spearman: pd.DataFrame, title: str, path: Path) -> None:
    """Spearman heatmap, wall time and peak memory of variant R1."""
    resources = summary[summary["ref"] == SPEARMAN_VARIANT].set_index("tool").loc[TOOLS]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), layout="constrained",
                             gridspec_kw={"width_ratios": [1.1, 1, 1]})
    draw_heatmap(fig, axes[0], spearman)
    draw_bars(axes[1], resources["wall_s"], "Mapper wall time", "s", 1)
    draw_bars(axes[2], resources["max_rss_mb"], "Peak memory", "MB", 0)
    fig.suptitle(f"{title}: host gene counts and mapper resources, variant {SPEARMAN_VARIANT}",
                 x=0.01, horizontalalignment="left", fontsize=10, color=INK_PRIMARY)
    add_caption(fig, ["One run per tool; wall time includes waiting on samtools sort for the aligners."])
    save(fig, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--resultdir", type=Path, required=True, help="results/benchmark/<dataset>")
    parser.add_argument("--title", required=True, help='figure title prefix, e.g. "Dataset A, SRR29668056"')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = pd.read_csv(args.resultdir / "benchmark_summary.tsv", sep="\t")
    summary["pct_end_clipped"] = 100 * summary["end_clipped_reads"] / summary["obelisk_alignments"]
    spearman = pd.read_csv(args.resultdir / f"host_spearman_{SPEARMAN_VARIANT}.tsv", sep="\t", index_col="tool")

    figdir = args.resultdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    plot_abundance(summary, args.title, figdir / "obelisk_abundance.png")
    plot_junction(summary, args.title, figdir / "obelisk_junction.png")
    plot_mapq_softclip(summary, args.title, figdir / "obelisk_mapq_softclip.png")
    plot_host_resources(summary, spearman, args.title, figdir / "host_spearman_resources.png")


if __name__ == "__main__":
    main()
