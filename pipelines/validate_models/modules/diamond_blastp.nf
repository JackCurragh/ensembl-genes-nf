// DIAMOND_BLASTP
// Translate predicted CDS to protein and DIAMOND blastp vs UniProt reviewed.
// Adds blast_pid and blast_cov attributes to the scored GFF3.
//
// DIAMOND is ~100× faster than BLASTP at equivalent sensitivity, making
// per-assembly protein validation feasible in minutes rather than hours.
//
// Inputs:
//   scored GFF3 from VALIDATE_MODELS (structural + splice scores already set)
//   genome FASTA (to extract CDS sequences)
//   protein DB FASTA (UniProt reviewed; converted to DIAMOND db on the fly)
//
// Outputs:
//   GFF3 with blast_pid + blast_cov + blast_desc attributes
//   TSV of per-protein DIAMOND top hits

process DIAMOND_BLASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::diamond=2.1.9 bioconda::samtools=1.18 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a88faf32df-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a88faf32df-0' }"

    input:
    tuple val(meta), path(scored_gff3)  // GFF3 from VALIDATE_MODELS
    path genome_fasta                    // softmasked genome (fai indexed)
    path protein_db_fasta               // UniProt reviewed FASTA

    output:
    tuple val(meta), path("*.validated.gff3"),   emit: gff3
    tuple val(meta), path("*.diamond_hits.tsv"), emit: hits
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def evalue = params.blast_evalue   ?: '1e-5'
    def min_id = params.blast_min_pid  ?: 30
    def min_cov= params.blast_min_cov  ?: 50
    def threads= task.cpus
    """
    # Build DIAMOND database from protein FASTA
    diamond makedb \\
        --in   ${protein_db_fasta} \\
        --db   uniprot.dmnd \\
        --threads ${threads}

    # Extract CDS sequences → translate to protein
    extract_cds_proteins.py \\
        --gff3   ${scored_gff3} \\
        --genome ${genome_fasta} \\
        --out    predicted_proteins.fa

    # DIAMOND blastp
    diamond blastp \\
        --db          uniprot.dmnd \\
        --query       predicted_proteins.fa \\
        --out         ${prefix}.diamond_hits.tsv \\
        --outfmt      6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle \\
        --evalue      ${evalue} \\
        --max-target-seqs 1 \\
        --threads     ${threads} \\
        ${args}

    # Annotate GFF3 with DIAMOND scores
    annotate_blast_scores.py \\
        --gff3      ${scored_gff3} \\
        --hits      ${prefix}.diamond_hits.tsv \\
        --min-pid   ${min_id} \\
        --min-cov   ${min_cov} \\
        --out       ${prefix}.validated.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | head -1 | sed 's/diamond version //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    cp ${scored_gff3} ${prefix}.validated.gff3
    printf 'qseqid\\tsseqid\\tpident\\tlength\\tqlen\\tslen\\tqcovhsp\\tevalue\\tbitscore\\tstitle\\n' \\
        > ${prefix}.diamond_hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: 2.1.9
        samtools: 1.18
        python: 3.11.0
    END_VERSIONS
    """
}
