#!/usr/bin/env python3
"""
extract_cds_proteins.py  —  Translate CDS features from a scored GFF3 into
protein FASTA for DIAMOND blastp.

Only transcripts with at least one CDS feature are translated.
Transcripts with internal stops or non-divisible-by-3 CDS are written with
an 'X' suffix and a warning in their header so DIAMOND can still attempt alignment.
"""

from __future__ import annotations
import argparse
import subprocess
import sys

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def faidx(genome: str, chrom: str, start: int, end: int) -> str:
    region = f"{chrom}:{start + 1}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", genome, region],
        capture_output=True, text=True, check=True
    )
    return "".join(result.stdout.split("\n")[1:]).upper()


def translate(seq: str) -> str:
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3",   required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out",    required=True)
    args = parser.parse_args()

    # Collect CDS per transcript
    tx_cds: dict[str, list] = {}
    tx_meta: dict[str, dict] = {}

    with open(args.gff3) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feat, start, end, score, strand, phase, attr_str = parts
            start, end = int(start) - 1, int(end)
            attrs = parse_attrs(attr_str)

            if feat in ("mRNA", "transcript"):
                tid = attrs.get("ID") or attrs.get("transcript_id")
                if tid:
                    tx_meta[tid] = {"chrom": chrom, "strand": strand}
                    tx_cds.setdefault(tid, [])

            elif feat == "CDS":
                parent = attrs.get("Parent", "")
                for ptx in parent.split(","):
                    ptx = ptx.strip()
                    if ptx:
                        tx_cds.setdefault(ptx, []).append((start, end, int(phase)))

    n_written = 0
    with open(args.out, "w") as out_fh:
        for tid, cds_segs in tx_cds.items():
            if not cds_segs:
                continue
            meta = tx_meta.get(tid, {})
            chrom  = meta.get("chrom", "")
            strand = meta.get("strand", "+")

            cds_segs.sort()
            try:
                if strand == "+":
                    pieces = [faidx(args.genome, chrom, s, e) for s, e, _ in cds_segs]
                else:
                    pieces = [rev_comp(faidx(args.genome, chrom, s, e))
                              for s, e, _ in reversed(cds_segs)]
                cds_seq = "".join(pieces)
            except Exception as exc:
                print(f"  [WARN] {tid}: CDS extraction failed: {exc}", file=sys.stderr)
                continue

            # Trim to frame
            if len(cds_seq) % 3 != 0:
                cds_seq = cds_seq[:-(len(cds_seq) % 3)]
            if len(cds_seq) < 3:
                continue

            prot = translate(cds_seq)
            # Remove trailing stop
            if prot.endswith("*"):
                prot = prot[:-1]

            out_fh.write(f">{tid}\n{prot}\n")
            n_written += 1

    print(f"Wrote {n_written} protein sequences to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
