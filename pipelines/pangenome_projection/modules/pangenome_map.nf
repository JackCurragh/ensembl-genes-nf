// PANGENOME_MAP
// Project gene models from a reference genome onto a target assembly using
// the Ensembl pangenome projection tool (minimap2 + cs-tag coordinate mapping).
//
// Tool: https://github.com/Ensembl/ensembl-genes/tree/feature/human_pangenome_mapping
//       pipelines/human_pangenome_projection/cli.py map
//
// The tool runs minimap2 internally (asm-mode PAF), parses cs-tags for exact
// coordinate projection of every exon/CDS/UTR, applies structural validation
// (splice sites, codon frame), and optional protein QC vs reference translations.
//
// Key output GFF3 attributes added vs chain-based projection:
//   mapped_from, mapping_identity, mapping_status, protein_identity
//
// Dependencies: minimap2, python>=3.11, pysam, biopython, edlib, click,
//               tqdm, rich, pandas, numpy, requests

process PANGENOME_MAP {
    tag "${meta.id}"
    label 'process_high'

    // conda spec covers all python deps + minimap2.
    // Container: use mulled or a custom image built with these deps.
    conda "bioconda::minimap2=2.26 bioconda::samtools=1.18 bioconda::pysam=0.22.0 conda-forge::python=3.11 conda-forge::biopython=1.82 conda-forge::edlib=1.3.9 conda-forge::click conda-forge::tqdm conda-forge::rich conda-forge::pandas conda-forge::numpy conda-forge::requests"

    input:
    tuple val(meta), path(source_fasta), path(source_gff3)   // reference genome + annotation
    path target_fasta                                          // target (new) assembly FASTA
    val  tool_dir                                              // path to cloned ensembl-genes repo

    output:
    tuple val(meta), path("*.pangenome_projected.gff3"), emit: gff3
    tuple val(meta), path("*.mapping_stats.json"),       emit: stats
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix         = task.ext.prefix ?: meta.id
    def args           = task.ext.args   ?: ''
    def preset         = params.pangenome_preset                   ?: 'none'
    def min_id         = params.pangenome_min_identity             ?: 0.95
    def min_block      = params.pangenome_min_block_length         ?: 10000
    def min_mapq       = params.pangenome_min_mapq                 ?: 10
    def pqc_min_id     = params.pangenome_protein_qc_min_identity  ?: 0.90
    def pqc_min_cov    = params.pangenome_protein_qc_min_coverage  ?: 0.90
    def threads        = task.cpus
    // tool_dir can be a pre-cloned path (preferred) or null → clone at runtime
    def setup_tool = tool_dir
        ? "PG_TOOL='${tool_dir}/pipelines/human_pangenome_projection'"
        : """
git clone --depth 1 --branch feature/human_pangenome_mapping \\
    https://github.com/Ensembl/ensembl-genes.git _ensembl_genes_pg
PG_TOOL="\${PWD}/_ensembl_genes_pg/pipelines/human_pangenome_projection"
"""
    """
    # ── Locate pangenome tool ─────────────────────────────────────────────────
    ${setup_tool}

    pip install -q pysam biopython edlib click tqdm rich pandas numpy requests 2>/dev/null || true
    export PYTHONPATH="\${PG_TOOL}:\${PYTHONPATH:-}"

    # ── Run projection ────────────────────────────────────────────────────────
    python "\${PG_TOOL}/cli.py" map \\
        --ref-fasta    ${source_fasta} \\
        --ref-gff      ${source_gff3} \\
        --target-fasta ${target_fasta} \\
        --output-gff   ${prefix}.pangenome_projected.gff3 \\
        --output-stats ${prefix}.mapping_stats.json \\
        --preset-profile ${preset} \\
        --min-identity   ${min_id} \\
        --min-block-length ${min_block} \\
        --min-mapq       ${min_mapq} \\
        --protein-qc-min-identity ${pqc_min_id} \\
        --protein-qc-min-coverage ${pqc_min_cov} \\
        --threads        ${threads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        python: \$(python3 --version | sed 's/Python //')
        ensembl_pangenome_tool: feature/human_pangenome_mapping
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\\n' > ${prefix}.pangenome_projected.gff3
    printf 'chr1\\tpangenome\\tgene\\t1000\\t5000\\t.\\t+\\t.\\t' >> ${prefix}.pangenome_projected.gff3
    printf 'ID=pg_gene_1;Name=ENSG001;biotype=projected_transcript;' >> ${prefix}.pangenome_projected.gff3
    printf 'mapping_identity=0.99;mapping_status=mapped\\n' >> ${prefix}.pangenome_projected.gff3

    echo '{"mapped": 100, "partial": 10, "unmapped": 5}' > ${prefix}.mapping_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.26
        python: 3.11.0
        ensembl_pangenome_tool: feature/human_pangenome_mapping
    END_VERSIONS
    """
}
