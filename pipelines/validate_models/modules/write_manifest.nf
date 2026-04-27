// WRITE_MANIFEST
// Collect all validated/scored GFF3 paths and write output_manifest.json.
// Unlike the consolidate manifest (single GFF3), validate_models may produce
// N scored GFF3s (one per annotation source), so this module accepts a list.

process WRITE_MANIFEST {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'quay.io/biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${outdir}/validate_models", mode: 'copy', overwrite: true

    input:
    val  outdir
    path gff3_files    // one or more validated GFF3 files staged in the work dir

    output:
    path 'output_manifest.json', emit: manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    // gff3_files may be a single path object or a list — normalise to space-sep
    def gff3_str = (gff3_files instanceof List ? gff3_files : [gff3_files])
        .collect { it.toString() }.join(' ')
    """
    write_manifest.py \\
        --pipeline validate_models \\
        --outdir   ${outdir}/validate_models \\
        --gff3s    ${gff3_str}
    """

    stub:
    """
    echo '{"pipeline":"validate_models","version":"1.0.0","completed_at":"stub","outputs":[]}' \
        > output_manifest.json
    """
}
