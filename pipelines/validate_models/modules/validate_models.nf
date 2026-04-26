// VALIDATE_MODELS
// Score gene models from any annotation source before consolidation.
//
// Checks applied to every transcript:
//   1. Structural sanity
//      - Splice sites: GT-AG / GC-AG / AT-AC vs non-canonical
//      - CDS: no in-frame stop codons, length divisible by 3
//      - ORF completeness: has start ATG + stop codon (if CDS annotated)
//   2. Splice junction support (only when --sj_tabs are provided)
//      - Fraction of predicted introns with ≥1 uniquely-mapped STAR read
//      - Median unique-read depth across supported junctions
//   3. Composite score (0–1) stored as GFF3 attribute
//      validation_score = 0.4 * structural_ok + 0.6 * splice_support
//      (if no RNA-seq junctions: validation_score = structural_ok)
//
// Outputs a scored GFF3 (same records, extra attributes added) and a TSV
// summary that the consolidation step can use to break layer-priority ties.
//
// Dependencies: samtools (sequence extraction), python ≥3.11

process VALIDATE_MODELS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.18 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0' }"

    input:
    tuple val(meta), path(gff3)         // GFF3 from any source sub-pipeline
    path genome_fasta                    // softmasked genome (must have .fai index)
    path sj_tabs                         // optional: list of STAR SJ.out.tab files (pass [] to skip)

    output:
    tuple val(meta), path("*.scored.gff3"), emit: gff3
    tuple val(meta), path("*.scores.tsv"),  emit: scores
    path "*.rejected.tsv",                  emit: rejected, optional: true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: meta.id
    def sj_arg        = sj_tabs ? "--sj-tabs ${(sj_tabs instanceof List ? sj_tabs : [sj_tabs]).join(' ')}" : ''
    def min_depth     = params.min_junction_depth        ?: 3
    def struct_pass   = params.structural_score_pass     ?: 0.7
    def struct_weight = params.validation_struct_weight  ?: 0.4
    def splice_weight = params.validation_splice_weight  ?: 0.6
    def rejected_arg  = "--rejected-tsv ${prefix}.rejected.tsv"
    """
    # Index genome if not already indexed
    if [ ! -f ${genome_fasta}.fai ]; then
        samtools faidx ${genome_fasta}
    fi

    validate_models.py \\
        --gff3             ${gff3} \\
        --genome           ${genome_fasta} \\
        --out-gff3         ${prefix}.scored.gff3 \\
        --out-tsv          ${prefix}.scores.tsv \\
        --source-label     ${prefix} \\
        --min-depth        ${min_depth} \\
        --structural-pass  ${struct_pass} \\
        --struct-weight    ${struct_weight} \\
        --splice-weight    ${splice_weight} \\
        ${rejected_arg} \\
        ${sj_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    cp ${gff3} ${prefix}.scored.gff3
    printf 'transcript_id\\tgene_id\\tsource\\tn_introns\\tcanonical_pct\\torf_complete\\tsplice_support_pct\\tmedian_depth\\tvalidation_score\\n' \\
        > ${prefix}.scores.tsv
    printf 'ENST00000000001\\tENSG00000000001\\t${prefix}\\t5\\t100.0\\ttrue\\t80.0\\t12.0\\t0.88\\n' \\
        >> ${prefix}.scores.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        samtools: 1.18
    END_VERSIONS
    """
}
