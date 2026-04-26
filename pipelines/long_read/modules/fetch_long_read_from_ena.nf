// FETCH_READS_FROM_ENA (long-read)
// Download long-read RNA-seq data from ENA by BioProject accession(s).
// Filters for OXFORD_NANOPORE and PACBIO_SMRT platforms with TRANSCRIPTOMIC
// library source. Produces a TSV sample sheet for the long_read pipeline.
//
// Accepts comma-separated BioProject accessions (e.g. PRJEB1234,PRJEB5678)
// so that runs from multiple projects can be combined in one pipeline run.
//
// Output TSV columns: sample_name, fastq_file, instrument_platform

process FETCH_READS_FROM_ENA {
    tag "${accession}"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::requests=2.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    val accession   // BioProject ID(s), comma-separated

    output:
    path "long_read_samples.tsv", emit: sample_sheet
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def max_runs  = params.max_longread_runs ?: 200
    """
    fetch_reads_from_ena.py \\
        --accession   "${accession}" \\
        --outdir      . \\
        --max-runs    ${max_runs} \\
        --platform    OXFORD_NANOPORE PACBIO_SMRT \\
        --output-tsv  long_read_samples.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """

    stub:
    """
    printf 'sample_name\\tfastq_file\\tinstrument_platform\\n' > long_read_samples.tsv
    printf 'spleen\\tSRR21742646_subreads.fastq.gz\\tPACBIO_SMRT\\n' >> long_read_samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        requests: 2.31.0
    END_VERSIONS
    """
}
