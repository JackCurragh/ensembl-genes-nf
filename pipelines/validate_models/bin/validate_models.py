#!/usr/bin/env python3
"""
validate_models.py  —  Score gene models before consolidation.

Checks applied to every transcript
───────────────────────────────────
1. Structural sanity
   a. Splice-site dinucleotides: canonical GT-AG / GC-AG / AT-AC vs non-canonical
   b. CDS completeness:
      • no in-frame stop codons
      • total CDS length divisible by 3
      • starts with ATG (if 5′-complete)
      • ends with stop codon (if 3′-complete)
2. Splice junction support (optional; requires --sj-tabs)
   • Fraction of predicted introns found in merged STAR SJ.out.tab
   • Median unique-read depth of supported junctions
3. Composite validation score (stored as GFF3 attribute)

Scoring formula
───────────────
  structural_score = 1.0 if all CDS + splice checks pass, else 0.0-0.8 depending on severity
  splice_support   = n_supported / n_introns  (1.0 if monoexonic or no junction data)
  validation_score = round(0.4 * structural_score + 0.6 * splice_support, 3)
                   = structural_score alone if no junction data supplied

Outputs
───────
  <prefix>.scored.gff3  : input GFF3 with validation attributes appended
  <prefix>.scores.tsv   : per-transcript summary (one row per transcript)

GFF3 attribute keys added
─────────────────────────
  validation_score      float 0-1
  structural_ok         true|false
  canonical_splice_pct  float 0-100
  orf_complete          true|false (false if no CDS features)
  n_introns             int
  splice_support_pct    float 0-100 (absent if no junction data)
  median_junc_depth     float (absent if no junction data)
"""

from __future__ import annotations

import argparse
import collections
import subprocess
import sys
from pathlib import Path

# ─── Genetic code (standard) ─────────────────────────────────────────────────
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

# Canonical splice site pairs (donor_dinuc, acceptor_dinuc) for both strands.
# We always read the genome in the 5′→3′ direction of the transcript.
CANONICAL   = {"GTAG", "GCAG", "ATAC"}   # merged donor+acceptor as 4-char key
NEAR_CANON  = {"GCAG"}                    # degraded but seen in real genomes

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


# ─── FASTA sequence fetcher (calls samtools faidx) ───────────────────────────
def faidx(genome: str, chrom: str, start: int, end: int) -> str:
    """
    Fetch sequence from indexed FASTA using samtools faidx.
    Coordinates are 0-based half-open [start, end).
    Returns upper-case sequence string.
    """
    region = f"{chrom}:{start + 1}-{end}"   # samtools is 1-based inclusive
    result = subprocess.run(
        ["samtools", "faidx", genome, region],
        capture_output=True, text=True, check=True
    )
    seq = "".join(result.stdout.split("\n")[1:]).upper()
    return seq


def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


# ─── GFF3 parser ──────────────────────────────────────────────────────────────
def parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def attrs_to_str(attrs: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in attrs.items())


def parse_gff3(path: str):
    """
    Parse a GFF3 file.  Returns:
      genes:       {gene_id → dict with children list}
      transcripts: {tx_id → dict(chrom, strand, exons, cds, gene_id, raw_attrs, raw_line)}
      raw_lines:   list of all original lines (for re-writing)
    """
    transcripts: dict = {}
    gene_of: dict[str, str] = {}   # transcript_id → gene_id
    raw_lines: list[str] = []

    with open(path) as fh:
        for line in fh:
            raw_lines.append(line)
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feat, start, end, score, strand, phase, attr_str = parts
            start, end = int(start) - 1, int(end)   # convert to 0-based half-open
            attrs = parse_attrs(attr_str)
            tid = attrs.get("ID") or attrs.get("transcript_id")
            parent = attrs.get("Parent", "")

            if feat in ("mRNA", "transcript"):
                transcripts[tid] = {
                    "chrom":  chrom,
                    "strand": strand,
                    "start":  start,
                    "end":    end,
                    "gene_id": parent,
                    "exons": [],
                    "cds":   [],
                    "raw_attrs": attrs,
                    "source": source,
                }
                gene_of[tid] = parent

            elif feat == "exon":
                for ptx in parent.split(","):
                    if ptx in transcripts:
                        transcripts[ptx]["exons"].append((start, end))

            elif feat == "CDS":
                for ptx in parent.split(","):
                    if ptx in transcripts:
                        transcripts[ptx]["cds"].append((start, end, int(phase)))

    # Sort exons / CDS
    for tx in transcripts.values():
        tx["exons"].sort()
        tx["cds"].sort()

    return transcripts, raw_lines


# ─── Structural checks ────────────────────────────────────────────────────────
def check_splice_sites(tx: dict, genome: str) -> tuple[int, int, float]:
    """
    Returns (n_introns, n_canonical, canonical_pct).
    Introns are derived from sorted exon coordinates.
    """
    exons = sorted(tx["exons"])
    if len(exons) < 2:
        return 0, 0, 100.0

    n_introns   = 0
    n_canonical = 0
    strand      = tx["strand"]
    chrom       = tx["chrom"]

    for i in range(len(exons) - 1):
        intron_start = exons[i][1]      # 0-based, end of exon i
        intron_end   = exons[i + 1][0]  # 0-based, start of exon i+1
        if intron_end <= intron_start:
            continue
        n_introns += 1
        try:
            donor    = faidx(genome, chrom, intron_start,     intron_start + 2)
            acceptor = faidx(genome, chrom, intron_end - 2,   intron_end)
            if strand == "-":
                donor, acceptor = rev_comp(acceptor), rev_comp(donor)
            key = donor + acceptor   # e.g. "GTAG"
            if key in CANONICAL:
                n_canonical += 1
        except Exception:
            pass   # network/index error — don't penalise

    pct = (100.0 * n_canonical / n_introns) if n_introns else 100.0
    return n_introns, n_canonical, pct


def get_intron_coords(tx: dict) -> list[tuple[str, int, int, str]]:
    """Return list of (chrom, start, end, strand) for each intron (0-based half-open)."""
    exons  = sorted(tx["exons"])
    chrom  = tx["chrom"]
    strand = tx["strand"]
    introns = []
    for i in range(len(exons) - 1):
        s = exons[i][1]
        e = exons[i + 1][0]
        if e > s:
            introns.append((chrom, s, e, strand))
    return introns


def check_cds(tx: dict, genome: str) -> tuple[bool, bool]:
    """
    Returns (orf_complete, no_internal_stops).
    orf_complete = has start ATG + ends with stop (if CDS present).
    """
    cds_segs = sorted(tx["cds"])
    if not cds_segs:
        return False, True   # no CDS → incompleteness, but not a stop-codon error

    strand = tx["strand"]
    chrom  = tx["chrom"]

    try:
        if strand == "+":
            pieces = [faidx(genome, chrom, s, e) for s, e, _ in cds_segs]
        else:
            pieces = [rev_comp(faidx(genome, chrom, s, e)) for s, e, _ in reversed(cds_segs)]
        cds_seq = "".join(pieces)
    except Exception:
        return False, True

    if len(cds_seq) % 3 != 0 or len(cds_seq) < 3:
        return False, True

    codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq), 3)]
    no_internal = all(c not in STOP_CODONS for c in codons[:-1])
    has_start   = codons[0] == START_CODON
    has_stop    = codons[-1] in STOP_CODONS
    orf_complete = has_start and has_stop

    return orf_complete, no_internal


def compute_structural_score(canonical_pct: float, orf_complete: bool, no_internal_stops: bool) -> float:
    """
    Return a structural score 0-1.
    Full credit (1.0) requires:
      • all splice sites canonical
      • ORF complete (start + stop)
      • no internal stops
    Penalties:
      • non-canonical splice site: -0.1 per 10 % below 100 %
      • missing start or stop: -0.15 each
      • internal stop: hard cap at 0.0
    """
    if not no_internal_stops:
        return 0.0

    score = 1.0
    # Penalise non-canonical splice sites
    if canonical_pct < 100.0:
        score -= 0.01 * (100.0 - canonical_pct)   # 1 % per percentage point below 100
    # Penalise incomplete ORF
    if not orf_complete:
        score -= 0.2
    return max(0.0, min(1.0, round(score, 3)))


# ─── Splice junction support ──────────────────────────────────────────────────
def load_sj_tabs(paths: list[str]) -> dict[tuple, int]:
    """
    Merge STAR SJ.out.tab files.
    Key: (chrom, intron_start_0based, intron_end_0based, strand_char)
    Value: total unique reads across all samples.
    STAR SJ.out.tab columns (1-based coords, inclusive):
      0  chrom
      1  start (1-based)
      2  end   (1-based, inclusive)
      3  strand (0=undefined, 1=+, 2=-)
      6  unique reads
    """
    junctions: dict[tuple, int] = collections.defaultdict(int)
    strand_map = {"0": ".", "1": "+", "2": "-"}

    for path in paths:
        try:
            with open(path) as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 7:
                        continue
                    chrom  = parts[0]
                    start  = int(parts[1]) - 1   # convert to 0-based
                    end    = int(parts[2])        # inclusive→half-open: keep as-is
                    strand = strand_map.get(parts[3], ".")
                    unique = int(parts[6])
                    junctions[(chrom, start, end, strand)] += unique
        except Exception as exc:
            print(f"  [WARN] Could not parse SJ.out.tab {path}: {exc}", file=sys.stderr)

    return junctions


def score_junction_support(introns: list, junctions: dict, min_depth: int) -> tuple[float, float]:
    """
    Returns (splice_support_pct, median_depth).
    Introns: list of (chrom, start, end, strand).
    """
    if not introns:
        return 100.0, 0.0

    depths = []
    n_supported = 0
    for chrom, start, end, strand in introns:
        # Try both the exact strand and '.' (undefined strand in SJ.out.tab)
        depth = junctions.get((chrom, start, end, strand), 0) + \
                junctions.get((chrom, start, end, "."), 0)
        depths.append(depth)
        if depth >= min_depth:
            n_supported += 1

    pct    = 100.0 * n_supported / len(introns)
    sorted_d = sorted(depths)
    n = len(sorted_d)
    if n % 2 == 0:
        median = (sorted_d[n // 2 - 1] + sorted_d[n // 2]) / 2
    else:
        median = float(sorted_d[n // 2])

    return round(pct, 1), round(median, 1)


# ─── Composite score ──────────────────────────────────────────────────────────
def composite_score(structural: float, splice_pct: float, has_junctions: bool) -> float:
    if not has_junctions:
        return round(structural, 3)
    splice = splice_pct / 100.0
    return round(0.4 * structural + 0.6 * splice, 3)


# ─── GFF3 writer with injected attributes ─────────────────────────────────────
def write_scored_gff3(raw_lines: list[str], scores: dict[str, dict], out_path: str) -> None:
    """
    Re-write the GFF3, appending validation attributes to mRNA/transcript lines.
    scores: {transcript_id → validation_attr_dict}
    """
    with open(out_path, "w") as fh:
        for line in raw_lines:
            if line.startswith("#") or not line.strip():
                fh.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fh.write(line)
                continue
            feat = parts[2]
            if feat in ("mRNA", "transcript"):
                attrs = parse_attrs(parts[8])
                tid   = attrs.get("ID") or attrs.get("transcript_id")
                if tid and tid in scores:
                    attrs.update(scores[tid])
                parts[8] = attrs_to_str(attrs)
                fh.write("\t".join(parts) + "\n")
            else:
                fh.write(line)


# ─── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3",         required=True,  help="Input GFF3 to validate")
    parser.add_argument("--genome",       required=True,  help="Softmasked genome FASTA (samtools faidx indexed)")
    parser.add_argument("--out-gff3",     required=True,  help="Output scored GFF3")
    parser.add_argument("--out-tsv",      required=True,  help="Output per-transcript TSV")
    parser.add_argument("--source-label", default="unknown", help="Label for the source sub-pipeline")
    parser.add_argument("--sj-tabs",      nargs="*", default=[], help="STAR SJ.out.tab files (0 or more)")
    parser.add_argument("--min-depth",    type=int, default=3,
                        help="Minimum unique-read depth to count a junction as 'supported' (default 3)")
    args = parser.parse_args()

    print(f"Validating models in {args.gff3}", file=sys.stderr)

    # ── Parse GFF3
    transcripts, raw_lines = parse_gff3(args.gff3)
    print(f"  {len(transcripts)} transcripts loaded", file=sys.stderr)

    # ── Load junction data
    has_junctions = bool(args.sj_tabs)
    junctions: dict = {}
    if has_junctions:
        print(f"  Loading {len(args.sj_tabs)} SJ.out.tab file(s)...", file=sys.stderr)
        junctions = load_sj_tabs(args.sj_tabs)
        print(f"  {len(junctions)} unique junctions in pool", file=sys.stderr)

    # ── Score each transcript
    score_attrs: dict[str, dict] = {}
    tsv_rows: list[dict] = []

    for tid, tx in transcripts.items():
        # Structural checks
        n_introns, n_canon, canon_pct = check_splice_sites(tx, args.genome)
        orf_complete, no_internal     = check_cds(tx, args.genome)
        struct_score = compute_structural_score(canon_pct, orf_complete, no_internal)
        structural_ok = (struct_score > 0.7)

        # Splice junction support
        splice_pct   = 100.0
        median_depth = 0.0
        if has_junctions and n_introns > 0:
            introns   = get_intron_coords(tx)
            splice_pct, median_depth = score_junction_support(introns, junctions, args.min_depth)

        # Composite
        val_score = composite_score(struct_score, splice_pct, has_junctions and n_introns > 0)

        attrs = {
            "validation_score":    str(val_score),
            "structural_ok":       "true" if structural_ok else "false",
            "canonical_splice_pct":str(round(canon_pct, 1)),
            "orf_complete":        "true" if orf_complete else "false",
            "n_introns":           str(n_introns),
        }
        if has_junctions and n_introns > 0:
            attrs["splice_support_pct"] = str(splice_pct)
            attrs["median_junc_depth"]  = str(median_depth)

        score_attrs[tid] = attrs

        tsv_rows.append({
            "transcript_id":      tid,
            "gene_id":            tx.get("gene_id", ""),
            "source":             args.source_label,
            "n_introns":          n_introns,
            "canonical_pct":      round(canon_pct, 1),
            "orf_complete":       orf_complete,
            "no_internal_stops":  no_internal,
            "splice_support_pct": round(splice_pct, 1),
            "median_depth":       round(median_depth, 1),
            "structural_score":   struct_score,
            "validation_score":   val_score,
        })

    # ── Write outputs
    write_scored_gff3(raw_lines, score_attrs, args.out_gff3)

    with open(args.out_tsv, "w") as tsv:
        header = ("transcript_id\tgene_id\tsource\tn_introns\tcanonical_pct"
                  "\torf_complete\tno_internal_stops\tsplice_support_pct"
                  "\tmedian_depth\tstructural_score\tvalidation_score\n")
        tsv.write(header)
        for row in tsv_rows:
            tsv.write(
                f"{row['transcript_id']}\t{row['gene_id']}\t{row['source']}\t"
                f"{row['n_introns']}\t{row['canonical_pct']}\t"
                f"{row['orf_complete']}\t{row['no_internal_stops']}\t"
                f"{row['splice_support_pct']}\t{row['median_depth']}\t"
                f"{row['structural_score']}\t{row['validation_score']}\n"
            )

    # ── Summary
    n_ok    = sum(1 for r in tsv_rows if r["structural_score"] > 0.7)
    n_total = len(tsv_rows)
    if has_junctions:
        mean_support = sum(r["splice_support_pct"] for r in tsv_rows) / n_total if n_total else 0
        print(f"  Structural OK    : {n_ok}/{n_total}", file=sys.stderr)
        print(f"  Mean junction support: {mean_support:.1f}%", file=sys.stderr)
    else:
        print(f"  Structural OK    : {n_ok}/{n_total}", file=sys.stderr)
    print(f"  → {args.out_gff3}", file=sys.stderr)
    print(f"  → {args.out_tsv}",  file=sys.stderr)


if __name__ == "__main__":
    main()
