#!/usr/bin/env python3
"""Build the sample sheets of both datasets from ENA run metadata.

Dataset A (PRJNA1129866, Zheludev et al. 2024): RNA-seq runs only. The design is
read from ``library_name`` (e.g. ``ObP1-A_rRNAdep`` = strain, replicate, prep).

Dataset B (PRJNA937727): S. sanguinis monoculture runs only. The design is read
from ``sample_title`` (e.g. ``S.sanguinis in Glucose, oxic conditions replicate 1``);
co-culture titles do not match and are left out.

Outputs written to ``--outdir``:
    samples_A.tsv    run, sample, strain, replicate, prep, read_count
    samples_B.tsv    run, sample, sugar, oxygen, replicate, read_count
"""

import argparse
import re
from pathlib import Path

import pandas as pd

A_PATTERN = r"^(?P<strain>ObP1|ObN1|CTRL)-(?P<replicate>[ABC])_(?P<prep>totalRNA|rRNAdep|RNaseR)$"
B_PATTERN = (
    r"^S\.sanguinis in (?P<sugar>glucose|sucrose), "
    r"(?P<oxygen>oxic|anoxic) conditions replicate (?P<replicate>\d)$"
)
B_EXPECTED_RUNS = 10


def samples_a(ena: pd.DataFrame) -> pd.DataFrame:
    """Return the dataset A sample sheet (RNA-seq runs only)."""
    rnaseq = ena[ena["library_strategy"] == "RNA-Seq"]
    design = rnaseq["library_name"].str.extract(A_PATTERN)
    sheet = pd.concat([rnaseq[["run_accession", "library_name"]], design, rnaseq[["read_count"]]], axis=1)
    sheet = sheet.rename(columns={"run_accession": "run", "library_name": "sample"})
    sheet["strain"] = sheet["strain"].replace({"CTRL": "media_control"})
    return sheet.sort_values(["prep", "strain", "replicate"])


def samples_b(ena: pd.DataFrame) -> pd.DataFrame:
    """Return the dataset B sample sheet (monoculture runs only)."""
    design = ena["sample_title"].str.extract(B_PATTERN, flags=re.IGNORECASE)
    sheet = pd.concat([ena[["run_accession"]], design, ena[["read_count"]]], axis=1).dropna()
    sheet["sugar"] = sheet["sugar"].str.lower()
    if len(sheet) != B_EXPECTED_RUNS:
        raise ValueError(f"Expected {B_EXPECTED_RUNS} monoculture runs, found {len(sheet)}")
    sheet = sheet.rename(columns={"run_accession": "run"})
    sheet.insert(1, "sample", sheet["sugar"] + "_" + sheet["oxygen"] + "_" + sheet["replicate"])
    return sheet.sort_values(["oxygen", "sugar", "replicate"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ena-a", type=Path, required=True, help="ENA run metadata for dataset A")
    parser.add_argument("--ena-b", type=Path, required=True, help="ENA run metadata for dataset B")
    parser.add_argument("--outdir", type=Path, required=True, help="output directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    samples_a(pd.read_csv(args.ena_a, sep="\t")).to_csv(args.outdir / "samples_A.tsv", sep="\t", index=False)
    samples_b(pd.read_csv(args.ena_b, sep="\t")).to_csv(args.outdir / "samples_B.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
