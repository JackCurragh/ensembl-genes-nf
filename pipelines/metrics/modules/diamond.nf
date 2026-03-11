process DIAMOND {
    label 'process_high'

    container "community.wave.seqera.io/library/diamond:2.1.24--61a5af76160d103f"

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/Diamond_alignments", mode: 'copy'

    input:
    tuple val(meta), path(ref_protein_faa)

    output:
    tuple val(meta), path("${meta.id}.dmnd"),           emit: db
    path "versions.yml",                                emit: versions

    script:
    """
    diamond blastp \
        --query proteins.faa \
        --db uniprot.dmnd \
        --threads 16 \
        --evalue 1e-5 \
        --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length qstart qend qlen qcovhsp sstart send slen evalue bitscore \
        --out hits.tsv

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


    script:
    """
    diamond blastp \
    --query proteins.faa \
    --db uniprot.dmnd \
    --threads 16 \
    --evalue 1e-5 \
    --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident length qstart qend qlen qcovhsp sstart send slen evalue bitscore \
    --out hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_c: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_C.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_c: 1.0.0
    END_VERSIONS
    """
}
