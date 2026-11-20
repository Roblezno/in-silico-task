#!/usr/bin/env python3
"""Compare the obelisk metrics of this work with the published values.

Dataset A is compared with Zheludev et al. (2024), whose values come as mean and
SD per library preparation (resources/zheludev2024_values.tsv). Dataset B is
compared with Maddamsetti & You (2025), whose values come per run
(resources/maddamsetti2025_values.tsv).

Outputs (in --outdir):
  obelisk_summary.tsv     one row per group, values formatted as "mean ± sd"
                          (sample SD), ready to go into the report
  b_vs_maddamsetti.tsv    one row per run of dataset B, numeric, with the ratio
                          between our and the published obelisk percentage

The published percentages are not exactly our metrics: they are defined on RNA
reads or on mapped reads, see the note column of the summary.
"""

import argparse
from pathlib import Path

import pandas as pd

SOURCE_A = "Zheludev et al. 2024"
SOURCE_B = "Maddamsetti & You 2025"
RANK_TEXT = "rank among 2271 genes by TPM"
NOTE_A = (
    "Published: % of RNA reads (k-mers, PolyBench with SmarterStrandedFilter), Fig. 4D and 5B; "
    f"non-rRNA % derived from the published means. This work: read pairs, HISAT2 + pysam; {RANK_TEXT}."
)
NOTE_B = (
    "Published: % of mapped reads (breseq), rank by unstranded kallisto TPM among 2271 genes, "
    f"Tables 1 and 2. This work: read pairs, HISAT2 + pysam; {RANK_TEXT}."
)

A_PREPS = ["totalRNA", "rRNAdep", "RNaseR"]
B_CONDITIONS = ["oxic", "anoxic"]
# Output name -> column of the obelisk metrics tables.
MEASURES = {
    "rrna_pct": "pct_rrna_mapped",
    "obelisk_pct_mapped": "pct_obelisk_mapped",
    "obelisk_pct_non_rrna": "pct_obelisk_non_rrna",
    "antisense_pct": "pct_antisense",
}
COMPARISON_COLUMNS = [
    "run", "sample", "oxygen", "sugar", "replicate", "pct_mapped", "pct_mapped_published",
    "obelisk_pct_mapped", "obelisk_pct_mapped_published", "obelisk_ratio",
    "obelisk_rank", "obelisk_rank_published",
]


def format_mean(mean: float) -> str:
    """Format a mean with 2 decimals below 1 and 1 decimal otherwise."""
    return f"{mean:.2f}" if mean < 1 else f"{mean:.1f}"


def format_mean_sd(mean: float, sd: float) -> str:
    """Format "mean ± sd" with the decimals set by the mean."""
    decimals = 2 if mean < 1 else 1
    return f"{mean:.{decimals}f} ± {sd:.{decimals}f}"


def format_series(values: pd.Series) -> str:
    """Format the mean and sample SD (ddof=1) of a series."""
    return format_mean_sd(values.mean(), values.std(ddof=1))


def format_rank(ranks: pd.Series) -> str:
    """Format ranks as "1" when all are equal, otherwise as "min-max"."""
    low, high = int(ranks.min()), int(ranks.max())
    return str(low) if low == high else f"{low}-{high}"


def zheludev_published(zheludev: pd.DataFrame, prep: str) -> dict[str, str]:
    """Return the published values of one preparation, formatted per output column."""
    values = zheludev[zheludev["prep"] == prep].set_index("metric")

    def formatted(metric: str) -> str:
        return format_mean_sd(values.loc[metric, "mean"], values.loc[metric, "sd"])

    obelisk, rrna = values.loc["pct_obelisk", "mean"], values.loc["pct_rrna", "mean"]
    non_rrna = 100 * obelisk / (100 - rrna)
    return {
        "rrna_pct": formatted("pct_rrna"),
        "obelisk_pct_mapped": formatted("pct_obelisk"),
        "obelisk_pct_non_rrna": f"{format_mean(non_rrna)} (derived)",
        "antisense_pct": formatted("pct_antisense"),
    }


def maddamsetti_published(group: pd.DataFrame) -> dict[str, str]:
    """Return the published values of the runs of one B group, formatted per output column."""
    return {
        "obelisk_pct_mapped": format_series(group["pct_obelisk_mapped_published"]),
        "obelisk_rank": format_rank(group["obelisk_rank_published"]),
    }


def summary_row(name: str, ours: pd.DataFrame, source: str, published: dict[str, str], note: str) -> dict:
    """Return one summary row: our metrics next to the published ones (empty if not reported)."""
    row = {"group": name, "n": len(ours), "published_source": source}
    for output, column in MEASURES.items():
        row[output] = format_series(ours[column])
        row[f"{output}_published"] = published.get(output, "")
    row["obelisk_rank"] = format_rank(ours["obelisk_rank"])
    row["obelisk_rank_published"] = published.get("obelisk_rank", "")
    row["note"] = note
    return row


def join_maddamsetti(metrics_b: pd.DataFrame, maddamsetti: pd.DataFrame) -> pd.DataFrame:
    """Add the published values to the B runs, stopping if a run is missing."""
    missing = metrics_b.loc[~metrics_b["run"].isin(maddamsetti["run"]), "run"]
    if not missing.empty:
        raise SystemExit(f"Error: runs of dataset B missing in the Maddamsetti file: {', '.join(missing)}")
    published = maddamsetti.drop(columns="source").rename(columns={
        "obelisk_rank": "obelisk_rank_published",
        "pct_mapped": "pct_mapped_published",
        "pct_obelisk_mapped": "pct_obelisk_mapped_published",
    })
    return metrics_b.merge(published, on="run", how="left", validate="one_to_one")


def build_summary(metrics_a: pd.DataFrame, b: pd.DataFrame, zheludev: pd.DataFrame) -> pd.DataFrame:
    """Return the summary table: the ObP1 groups of A, then the oxygen groups of B."""
    rows = []
    for prep in A_PREPS:
        ours = metrics_a[(metrics_a["strain"] == "ObP1") & (metrics_a["prep"] == prep)]
        rows.append(summary_row(f"A ObP1 {prep}", ours, SOURCE_A, zheludev_published(zheludev, prep), NOTE_A))
    for oxygen in B_CONDITIONS:
        ours = b[b["oxygen"] == oxygen]
        rows.append(summary_row(f"B {oxygen}", ours, SOURCE_B, maddamsetti_published(ours), NOTE_B))
    return pd.DataFrame(rows)


def build_comparison(b: pd.DataFrame) -> pd.DataFrame:
    """Return the numeric per-run comparison of dataset B with Maddamsetti & You."""
    table = b.assign(
        obelisk_ratio=(b["pct_obelisk_mapped"] / b["pct_obelisk_mapped_published"]).round(2),
        obelisk_rank=b["obelisk_rank"].astype(int),
        obelisk_rank_published=b["obelisk_rank_published"].astype(int),
    )
    table = table.rename(columns={
        "pct_obelisk_mapped": "obelisk_pct_mapped",
        "pct_obelisk_mapped_published": "obelisk_pct_mapped_published",
    })
    return table[COMPARISON_COLUMNS]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics-a", type=Path, required=True, help="obelisk metrics of dataset A")
    parser.add_argument("--metrics-b", type=Path, required=True, help="obelisk metrics of dataset B")
    parser.add_argument("--zheludev", type=Path, required=True, help="published values, Zheludev et al. 2024")
    parser.add_argument("--maddamsetti", type=Path, required=True, help="published values, Maddamsetti & You 2025")
    parser.add_argument("--outdir", type=Path, required=True, help="output directory (created)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metrics_a = pd.read_csv(args.metrics_a, sep="\t")
    metrics_b = pd.read_csv(args.metrics_b, sep="\t")
    zheludev = pd.read_csv(args.zheludev, sep="\t")
    maddamsetti = pd.read_csv(args.maddamsetti, sep="\t")

    b = join_maddamsetti(metrics_b, maddamsetti)
    summary = build_summary(metrics_a, b, zheludev)
    comparison = build_comparison(b)

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.outdir / "obelisk_summary.tsv", sep="\t", index=False)
    comparison.to_csv(args.outdir / "b_vs_maddamsetti.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    ratio = comparison["obelisk_ratio"]
    print(f"\nobelisk_ratio (ours / published), min-max: {ratio.min():.2f}-{ratio.max():.2f}")


if __name__ == "__main__":
    main()
