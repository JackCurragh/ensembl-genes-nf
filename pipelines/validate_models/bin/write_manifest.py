#!/usr/bin/env python3
"""
write_manifest.py — Write output_manifest.json for the validate_models pipeline.

Unlike the consolidate manifest (single GFF3), validate_models emits one
scored GFF3 per annotation source, so this script accepts N --gff3s arguments.
"""

import argparse
import json
import os
from datetime import datetime, timezone


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pipeline', default='validate_models')
    ap.add_argument('--version',  default='1.0.0')
    ap.add_argument('--outdir',   required=True,
                    help='Publish directory (validate_models sub-dir of stage outdir)')
    ap.add_argument('--gff3s',    required=True, nargs='+',
                    help='Staged filenames of the validated GFF3 files')
    args = ap.parse_args()

    manifest = {
        'pipeline':     args.pipeline,
        'version':      args.version,
        'completed_at': datetime.now(timezone.utc).isoformat(),
        'outputs': [
            {'type': 'gff3', 'path': os.path.join(args.outdir, gff3)}
            for gff3 in args.gff3s
        ],
    }

    with open('output_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
