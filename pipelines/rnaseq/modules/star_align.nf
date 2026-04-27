// STAR_ALIGN
// Align paired-end (or single-end) RNA-seq reads with STAR.
// Produces sorted BAM + SJ.out.tab splice junction file.

process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::star=2.7.11b bioconda::samtools=1.21"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"

    publishDir "${params.outdir}/star_align/${meta.id}", mode: 'copy', pattern: "*.{bam,bai,tab}"

    input:
    tuple val(meta), path(reads)   // reads: [ R1.fq.gz ] or [ R1.fq.gz, R2.fq.gz ]
    path  index_dir               // STAR genome index directory

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"),    emit: bam
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("*.SJ.out.tab"),                        emit: sj
    tuple val(meta), path("*.Log.final.out"),                     emit: log
    path  "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: meta.id
    def args      = task.ext.args   ?: ''
    def reads_arg = reads instanceof List ? reads.join(' ') : reads
    def pe_flag   = (reads instanceof List && reads.size() > 1) ? '' : '--readFilesIn'
    """
    # Decompress reads if gzipped — bgzip reads both gzip and bgzip format.
    # This avoids --readFilesCommand quoting issues across Linux/macOS conda envs.
    decomp_reads=""
    for f in ${reads_arg}; do
        if [[ "\$f" == *.gz ]]; then
            out="\${f%.gz}"
            bgzip -d -c "\$f" > "\$out"
            decomp_reads="\$decomp_reads \$out"
        else
            decomp_reads="\$decomp_reads \$f"
        fi
    done

    STAR \\
        --runThreadN     ${task.cpus} \\
        --genomeDir      ${index_dir} \\
        --readFilesIn    \${decomp_reads} \\
        --outSAMtype     BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        --outSAMattributes NH HI AS NM \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNmax 10 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFileNamePrefix ${prefix}. \\
        --runRNGseed       0 \\
        ${args}

    samtools index ${prefix}.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam.bai
    touch ${prefix}.SJ.out.tab
    printf 'STAR   Mapping speed, Million of reads per hour\t0\\n' > ${prefix}.Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11b
        samtools: 1.21
    END_VERSIONS
    """
}
