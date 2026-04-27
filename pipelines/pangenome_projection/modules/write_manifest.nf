// WRITE_MANIFEST — write output_manifest.json for pangenome_projection pipeline.
// Single GFF3 output (one projection per target assembly).

process WRITE_MANIFEST {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'quay.io/biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${outdir}/pangenome_projection", mode: 'copy', overwrite: true

    input:
    val  outdir
    path gff3

    output:
    path 'output_manifest.json', emit: manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - <<'PYEOF'
import json, os
from datetime import datetime, timezone

manifest = {
    'pipeline':     'pangenome_projection',
    'version':      '1.0.0',
    'completed_at': datetime.now(timezone.utc).isoformat(),
    'outputs': [
        {'type': 'gff3', 'path': os.path.join('${outdir}/pangenome_projection', '${gff3}')}
    ],
}
with open('output_manifest.json', 'w') as fh:
    json.dump(manifest, fh, indent=2)
print(json.dumps(manifest, indent=2))
PYEOF
    """

    stub:
    """
    echo '{"pipeline":"pangenome_projection","version":"1.0.0","completed_at":"stub","outputs":[]}' \
        > output_manifest.json
    """
}
