process PREP_AA_FASTA {
    label 'process_medium'

    container "community.wave.seqera.io/library/gffread:0.12.7--33b95f1cfcc0e572"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/tool_a",
        mode: 'copy',
        pattern: "*_A.txt"

    input:
    tuple val(meta), path(genome_fa)
    tuple val(meta2), path(annotation_gff3)

    output:
    tuple val(meta), path("${meta.id}.faa"),            emit: aa_fasta
    path "versions.yml",                                emit: versions

    script:
    """
    gffread ${annotation_gff3} -g ${genome_fa} -x cds.fa -y ${meta.id}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: $(gffread --version 2>&1 | sed 's/^gffread v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: $(gffread --version 2>&1 | sed 's/^gffread v//')
    END_VERSIONS
    """
}
