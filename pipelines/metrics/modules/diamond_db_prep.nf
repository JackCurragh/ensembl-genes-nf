process PREP_DIAMOND_DB {
    label 'process_medium'

    container "community.wave.seqera.io/library/diamond:2.1.24--61a5af76160d103f"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/Diamond_db", mode: 'copy'

    input:
    tuple val(meta), path(ref_protein_faa)

    output:
    tuple val(meta), path("${meta.id}.dmnd"),           emit: db
    path "versions.yml",                                emit: versions

    script:
    """
    diamond makedb \
        --in ${ref_protein_faa} \
        --db ${meta.id}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: $(diamond --version | sed 's/diamond version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: $(diamond --version | sed 's/diamond version //')
    END_VERSIONS
    """
}
