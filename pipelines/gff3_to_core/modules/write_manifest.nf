process WRITE_MANIFEST {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${outdir}", mode: 'copy', overwrite: true

    input:
    val  outdir
    path load_stats
    path stable_id_stats

    output:
    path 'output_manifest.json', emit: manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    write_manifest.py \\
        --pipeline  gff3_to_core \\
        --version   1.0.0 \\
        --outdir    ${outdir} \\
        --stats-json ${load_stats} \\
        --dbname    ${params.db_name}
    """

    stub:
    """
    echo '{"pipeline":"gff3_to_core","version":"1.0.0","completed_at":"stub","outputs":[{"type":"core_db","path":"stub","meta":{"dbname":"stub_db"}}]}' \\
        > output_manifest.json
    """
}
