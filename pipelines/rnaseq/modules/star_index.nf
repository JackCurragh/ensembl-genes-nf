// STAR_INDEX
// Build a STAR genome index from an unmasked genome FASTA.
// The index directory is published to outdir/star_index/.

process STAR_INDEX {
    label 'process_high'

    conda "bioconda::star=2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'quay.io/biocontainers/star:2.7.11b--h43eeafb_0' }"

    // Building a STAR index takes 30–60 min and 30+ GB RAM. Cache it so any
    // re-run of the RNA-seq pipeline (or a second RNA-seq experiment against
    // the same assembly) skips index generation entirely.
    // publishDir is intentionally omitted — storeDir IS the permanent location.
    storeDir "${params.outdir}/store/star_index"

    input:
    path genome_fasta
    path gtf           // optional annotation; pass assets/no_annotation.gtf (empty) to skip

    output:
    path "star_index/", emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    // Empty sentinel file (assets/no_annotation.gtf, size 0) means: skip junction DB.
    def gtf_arg  = (gtf.size() > 0) ? "--sjdbGTFfile ${gtf} --sjdbOverhang ${params.sjdb_overhang}" : ''
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        ${gtf_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/Genome star_index/SA star_index/SAindex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11b
    END_VERSIONS
    """
}
