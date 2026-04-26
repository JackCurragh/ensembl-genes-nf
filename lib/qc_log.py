"""
qc_log.py — Standardised rejection logging for ensembl-genes-nf filtering steps.

Every step that excludes data should call QCLog.reject() for each excluded item,
then QCLog.write() at the end.  The master pipeline collects all rejection TSVs
into a single QC report so you can audit why any run, transcript, or model was
dropped, and at what threshold.

TSV output columns
──────────────────
stage            Name of the pipeline stage writing this log
item_type        Type of thing being rejected: run | transcript | model | read
item_id          Identifier of the rejected item (run_accession, transcript_id, etc.)
reason           Short code for why it was rejected (see REASONS below)
threshold_name   Name of the parameter that set the cutoff
threshold_value  Value the threshold was set to
observed_value   The actual value of the rejected item (may be '' if not applicable)

Reason codes (extend as needed)
────────────────────────────────
low_mapping_rate        — STAR uniquely_mapped_reads_percentage below threshold
wrong_platform          — instrument_platform not in expected set
low_coverage            — alignment or assembly coverage below threshold
low_pid                 — protein identity below threshold
non_canonical_splice    — canonical_splice_pct below threshold
internal_stop           — in-frame stop codon found in CDS
incomplete_orf          — missing start or stop codon
low_splice_support      — splice_support_pct below threshold
low_stringtie_coverage  — StringTie cov attribute below threshold
short_transcript        — transcript length below minimum
few_exons               — exon count below minimum

Usage
─────
    from qc_log import QCLog

    log = QCLog("rnaseq_qc", output_path="rnaseq_qc.rejected.tsv")
    log.reject("run", "SRR123456", "low_mapping_rate",
               threshold_name="min_mapping_pct", threshold_value=60.0,
               observed_value=45.2)
    log.write()
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


COLUMNS = [
    "stage",
    "item_type",
    "item_id",
    "reason",
    "threshold_name",
    "threshold_value",
    "observed_value",
]


@dataclass
class Rejection:
    stage:           str
    item_type:       str
    item_id:         str
    reason:          str
    threshold_name:  str
    threshold_value: str
    observed_value:  str = ""


class QCLog:
    """Collect rejection events for one filtering step and write them to a TSV."""

    def __init__(self, stage: str, output_path: Optional[str] = None):
        self.stage       = stage
        self.output_path = output_path
        self._records: list[Rejection] = []

    def reject(
        self,
        item_type: str,
        item_id: str,
        reason: str,
        threshold_name: str,
        threshold_value,
        observed_value="",
    ) -> None:
        self._records.append(
            Rejection(
                stage=self.stage,
                item_type=item_type,
                item_id=str(item_id),
                reason=reason,
                threshold_name=threshold_name,
                threshold_value=str(threshold_value),
                observed_value=str(observed_value),
            )
        )

    def n_rejected(self) -> int:
        return len(self._records)

    def write(self, path: Optional[str] = None) -> None:
        out = path or self.output_path
        if not out:
            return
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
            writer.writeheader()
            for r in self._records:
                writer.writerow({
                    "stage":           r.stage,
                    "item_type":       r.item_type,
                    "item_id":         r.item_id,
                    "reason":          r.reason,
                    "threshold_name":  r.threshold_name,
                    "threshold_value": r.threshold_value,
                    "observed_value":  r.observed_value,
                })
        print(
            f"  QC log: {len(self._records)} rejected items → {out}",
            file=sys.stderr,
        )

    def print_summary(self) -> None:
        """Print a grouped summary to stderr."""
        if not self._records:
            return
        by_reason: dict[str, int] = {}
        for r in self._records:
            by_reason[r.reason] = by_reason.get(r.reason, 0) + 1
        print(f"\n  [{self.stage}] Rejection summary:", file=sys.stderr)
        for reason, count in sorted(by_reason.items(), key=lambda x: -x[1]):
            print(f"    {count:>6}  {reason}", file=sys.stderr)
