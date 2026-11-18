#!/usr/bin/env python3
"""Build the reference files for the host genome plus Obelisk-S.s.

Inputs are the host GenBank record (sequence and annotation) and the obelisk
FASTA. The obelisk is circular, so the first ``--extension`` nt are appended to
its end; reads that cross the circle junction by up to that length can then
align in one piece. The obelisk is annotated twice over the same contig: on
the + strand, the strand of its Oblin-1 ORF (sense), and on the - strand
(antisense). Stranded counting (featureCounts, kallisto) then quantifies each
strand separately; the antisense strand is the most abundant in Zheludev et
al. (2024).

Outputs written to ``--outdir``:
    sk36_obelisk.fa   host genome + extended obelisk (mapping, IGV)
    annotation.saf    host genes + obelisk sense and antisense (featureCounts)
    host_genes.bed    BED12 with host genes only (RSeQC infer_experiment.py)
    transcripts.fa    one sequence per gene, on its own strand, + obelisk
                      sense and antisense (kallisto)
    genes.tsv         gene table: ID, biotype, name, product, coordinates, length
"""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

OBELISK_BIOTYPE = "obelisk"


@dataclass
class Gene:
    """One annotated feature. Coordinates are 0-based, half-open."""

    gene_id: str
    biotype: str
    name: str
    product: str
    chrom: str
    start: int
    end: int
    strand: str


def first_qualifier(feature, key: str) -> str:
    """Return the first value of a GenBank qualifier, or an empty string."""
    return feature.qualifiers.get(key, [""])[0]


def collect_host_genes(record: SeqRecord) -> list[Gene]:
    """Return one Gene per ``gene`` feature, typed by its child feature.

    The child (CDS, rRNA, tRNA, misc_RNA) shares the gene's locus_tag and gives
    the biotype and product. Genes without a child are pseudogenes.
    """
    children = {
        first_qualifier(f, "locus_tag"): f
        for f in record.features
        if f.type not in ("source", "gene")
    }
    genes = []
    for feature in record.features:
        if feature.type != "gene":
            continue
        locus_tag = first_qualifier(feature, "locus_tag")
        child = children.get(locus_tag)
        genes.append(
            Gene(
                gene_id=locus_tag,
                biotype=child.type if child else "pseudogene",
                name=first_qualifier(feature, "gene"),
                product=first_qualifier(child, "product") if child else "",
                chrom=record.id,
                start=int(feature.location.start),
                end=int(feature.location.end),
                strand="+" if feature.location.strand == 1 else "-",
            )
        )
    return genes


def extend_circular(record: SeqRecord, extension: int) -> SeqRecord:
    """Append the first ``extension`` nt of a circular sequence to its end."""
    extended = record[:] + record[:extension]
    extended.id = record.id
    extended.description = f"{record.description} | extended +{extension} nt"
    return extended


def obelisk_genes(obelisk: SeqRecord) -> list[Gene]:
    """Annotate the whole (extended) obelisk contig as a sense (+) and an antisense (-) gene."""
    return [
        Gene(
            gene_id=obelisk.id,
            biotype=OBELISK_BIOTYPE,
            name="Obelisk-S.s",
            product="viroid-like circular RNA, sense strand (Oblin-1)",
            chrom=obelisk.id,
            start=0,
            end=len(obelisk),
            strand="+",
        ),
        Gene(
            gene_id=f"{obelisk.id}_antisense",
            biotype=OBELISK_BIOTYPE,
            name="Obelisk-S.s antisense",
            product="viroid-like circular RNA, antisense strand",
            chrom=obelisk.id,
            start=0,
            end=len(obelisk),
            strand="-",
        ),
    ]


def write_saf(genes: list[Gene], path: Path) -> None:
    """Write a featureCounts SAF file (1-based, inclusive coordinates)."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
        for gene in genes:
            writer.writerow([gene.gene_id, gene.chrom, gene.start + 1, gene.end, gene.strand])


def write_bed12(genes: list[Gene], path: Path) -> None:
    """Write a single-block BED12 file (bacterial genes have no introns)."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for gene in genes:
            length = gene.end - gene.start
            writer.writerow([
                gene.chrom, gene.start, gene.end, gene.gene_id, 0, gene.strand,
                gene.start, gene.end, 0, 1, f"{length},", "0,",
            ])


def write_transcripts(genes: list[Gene], sequences: dict[str, Seq], path: Path) -> None:
    """Write each gene's sequence on its own strand, named by gene ID."""
    records = []
    for gene in genes:
        segment = sequences[gene.chrom][gene.start:gene.end]
        if gene.strand == "-":
            segment = segment.reverse_complement()
        records.append(SeqRecord(segment, id=gene.gene_id, description=gene.biotype))
    SeqIO.write(records, path, "fasta")


def write_gene_table(genes: list[Gene], path: Path) -> None:
    """Write the gene table with 1-based, inclusive coordinates."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "biotype", "name", "product", "chrom", "start", "end", "strand", "length"])
        for gene in genes:
            writer.writerow([
                gene.gene_id, gene.biotype, gene.name, gene.product, gene.chrom,
                gene.start + 1, gene.end, gene.strand, gene.end - gene.start,
            ])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genbank", type=Path, required=True, help="host GenBank record")
    parser.add_argument("--obelisk", type=Path, required=True, help="obelisk FASTA (one sequence)")
    parser.add_argument("--extension", type=int, default=30, help="nt appended to the obelisk end (default: 30)")
    parser.add_argument("--outdir", type=Path, required=True, help="output directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    host = SeqIO.read(args.genbank, "genbank")
    obelisk = extend_circular(SeqIO.read(args.obelisk, "fasta"), args.extension)
    host_genes = collect_host_genes(host)
    all_genes = host_genes + obelisk_genes(obelisk)
    sequences = {host.id: host.seq, obelisk.id: obelisk.seq}

    SeqIO.write([host, obelisk], args.outdir / "sk36_obelisk.fa", "fasta")
    write_saf(all_genes, args.outdir / "annotation.saf")
    write_bed12(host_genes, args.outdir / "host_genes.bed")
    write_transcripts(all_genes, sequences, args.outdir / "transcripts.fa")
    write_gene_table(all_genes, args.outdir / "genes.tsv")


if __name__ == "__main__":
    main()
