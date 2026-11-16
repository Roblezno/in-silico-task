#!/usr/bin/env python3
"""Extract one obelisk sequence from the Zheludev et al. (2024) supplementary table.

The spreadsheet (``mmc2.xlsx``) lists every obelisk with its ID in ``Obelisk_ID``
and its sequence in ``Obelisk_seq``. This script runs once: its output FASTA is
versioned in ``resources/`` so the analysis does not depend on the spreadsheet.

Usage:
    python scripts/extract_obelisk.py --xlsx mmc2.xlsx \
        --obelisk-id Obelisk_000003_000001_000001 --name Obelisk_000003 \
        --output resources/obelisk_000003.fasta
"""

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SOURCE = "Zheludev et al. (2024) Cell 187:6521, supplementary table mmc2.xlsx"


def read_obelisk_sequence(xlsx: Path, obelisk_id: str) -> str:
    """Return the sequence of ``obelisk_id`` from the supplementary table."""
    table = pd.read_excel(xlsx, index_col="Obelisk_ID")
    return str(table.at[obelisk_id, "Obelisk_seq"]).upper()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xlsx", type=Path, required=True, help="supplementary table (mmc2.xlsx)")
    parser.add_argument("--obelisk-id", required=True, help="row ID in the Obelisk_ID column")
    parser.add_argument("--name", required=True, help="sequence name in the output FASTA")
    parser.add_argument("--output", type=Path, required=True, help="output FASTA")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sequence = read_obelisk_sequence(args.xlsx, args.obelisk_id)
    record = SeqRecord(
        Seq(sequence),
        id=args.name,
        description=f"{args.obelisk_id} | {SOURCE}",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(record, args.output, "fasta")


if __name__ == "__main__":
    main()
