#!/usr/bin/env python3
"""
annotate_blast_scores.py  —  Merge DIAMOND blastp results into a scored GFF3.

Reads the DIAMOND tabular output (outfmt 6 with:
  qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle)
and appends blast_pid, blast_cov, blast_desc, blast_validated attributes to
each mRNA/transcript line.

blast_validated = true if pident ≥ --min-pid AND qcovhsp ≥ --min-cov.
"""

from __future__ import annotations
import argparse
import sys


def parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def attrs_to_str(attrs: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in attrs.items())


def load_diamond_hits(path: str, min_pid: float, min_cov: float) -> dict[str, dict]:
    """Return {query_id → {pid, cov, desc, validated}}. Best hit only (first line)."""
    hits: dict[str, dict] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            qseqid, sseqid, pident, length, qlen, slen, qcov, evalue, bitscore, stitle = parts[:10]
            if qseqid in hits:
                continue   # keep best hit only
            pid = float(pident)
            cov = float(qcov)
            # Clean up the subject title (remove OS= etc for brevity)
            desc = stitle.split(" OS=")[0].strip()
            hits[qseqid] = {
                "pid":       round(pid, 1),
                "cov":       round(cov, 1),
                "desc":      desc.replace(";", ","),   # GFF3-safe
                "validated": pid >= min_pid and cov >= min_cov,
            }
    return hits


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3",    required=True)
    parser.add_argument("--hits",    required=True, help="DIAMOND tabular output")
    parser.add_argument("--out",     required=True)
    parser.add_argument("--min-pid", type=float, default=30.0)
    parser.add_argument("--min-cov", type=float, default=50.0)
    args = parser.parse_args()

    hits = load_diamond_hits(args.hits, args.min_pid, args.min_cov)
    print(f"Loaded {len(hits)} DIAMOND hits", file=sys.stderr)

    n_validated = 0
    with open(args.gff3) as fh, open(args.out, "w") as out_fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                out_fh.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                out_fh.write(line)
                continue
            feat = parts[2]
            if feat not in ("mRNA", "transcript"):
                out_fh.write(line)
                continue

            attrs = parse_attrs(parts[8])
            tid   = attrs.get("ID") or attrs.get("transcript_id")
            if tid and tid in hits:
                h = hits[tid]
                attrs["blast_pid"]       = str(h["pid"])
                attrs["blast_cov"]       = str(h["cov"])
                attrs["blast_desc"]      = h["desc"]
                attrs["blast_validated"] = "true" if h["validated"] else "false"

                # Combine blast into validation_score (if structural score present)
                if "validation_score" in attrs and h["validated"]:
                    old_score = float(attrs["validation_score"])
                    # Boost by up to 0.1 for blast support, capped at 1.0
                    boost = 0.1 * (h["pid"] / 100.0)
                    attrs["validation_score"] = str(round(min(1.0, old_score + boost), 3))
                    n_validated += 1

            parts[8] = attrs_to_str(attrs)
            out_fh.write("\t".join(parts) + "\n")

    print(f"  blast_validated=true for {n_validated} transcripts", file=sys.stderr)


if __name__ == "__main__":
    main()
