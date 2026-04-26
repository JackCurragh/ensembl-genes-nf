#!/usr/bin/env nextflow
/*
========================================================================================
    PANGENOME_PROJECTION PIPELINE
========================================================================================
    Projects gene models from a reference genome onto a target assembly using
    the Ensembl pangenome projection tool.

    Approach: minimap2 (asm-mode PAF) + exact cs-tag coordinate projection.
    Unlike the chain-based PROJECTION pipeline this does not require a pre-built
    lastz/liftOver chain — minimap2 aligns on the fly.

    Runs ALONGSIDE the chain-based projection for side-by-side comparison.
    Both outputs flow into VALIDATE_MODELS and CONSOLIDATE; the validation
    scores (structural + protein QC) will reveal which approach produces
    higher-quality models for each assembly.

    Inputs:
      --source_fasta   Reference genome FASTA (e.g. GRCh38)
      --source_gff3    Reference annotation GFF3 (e.g. Ensembl 110 human)
      --target_fasta   Target assembly FASTA (softmasked)
      --outdir         Output directory

    Optional:
      --pangenome_tool_dir  Path to a pre-cloned ensembl-genes repo at
                            feature/human_pangenome_mapping. If absent, the
                            process will git-clone at runtime (needs outbound
                            internet on compute nodes).

    Run locally (stub):
      nextflow run pipelines/pangenome_projection -profile local -stub \
        --source_fasta  grch38.fa \
        --source_gff3   grch38.gff3 \
        --target_fasta  target.fa \
        --outdir        results/pangenome_projection
========================================================================================
*/

nextflow.enable.dsl = 2

include { PANGENOME_MAP   } from './modules/pangenome_map.nf'
include { WRITE_MANIFEST  } from './modules/write_manifest.nf'

params.source_fasta        = null
params.source_gff3         = null
params.target_fasta        = null
params.pangenome_tool_dir  = null
params.outdir              = null

def validate_params() {
    def errors = []
    if (!params.source_fasta)  errors << "  --source_fasta is required"
    if (!params.source_gff3)   errors << "  --source_gff3 is required"
    if (!params.target_fasta)  errors << "  --target_fasta is required"
    if (!params.outdir)        errors << "  --outdir is required"
    if (!params.pangenome_tool_dir) {
        log.warn "  --pangenome_tool_dir not set; will clone ensembl-genes from GitHub at runtime"
        log.warn "  This requires outbound internet access on compute nodes."
        log.warn "  To avoid this: git clone --depth 1 --branch feature/human_pangenome_mapping"
        log.warn "  https://github.com/Ensembl/ensembl-genes.git /path/to/ensembl-genes"
        log.warn "  then pass --pangenome_tool_dir /path/to/ensembl-genes"
    }
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    // Assembly label: derive from target FASTA filename (strip path + extensions)
    def asm_label = file(params.target_fasta).name.replaceAll(/\.(softmasked\.fa|fa|fasta)(\.gz)?$/, '')
    ch_meta = Channel.value([id: asm_label])

    // Source reference: pass as a tuple so PANGENOME_MAP sees (meta, fasta, gff3)
    ch_source = ch_meta.map { meta ->
        tuple(meta,
              file(params.source_fasta, checkIfExists: true),
              file(params.source_gff3,  checkIfExists: true))
    }

    ch_target = file(params.target_fasta, checkIfExists: true)
    def tool_dir = params.pangenome_tool_dir ?: ''

    // ── Map ───────────────────────────────────────────────────────────────────
    PANGENOME_MAP(ch_source, ch_target, tool_dir)

    // ── Manifest ──────────────────────────────────────────────────────────────
    WRITE_MANIFEST(
        params.outdir,
        PANGENOME_MAP.out.gff3.map { meta, gff3 -> gff3 }
    )

    // Software versions
    Channel.empty()
        .mix(PANGENOME_MAP.out.versions)
        .collectFile(
            name:     'software_versions.tsv',
            newLine:  true,
            storeDir: "${params.outdir}/pipeline_info"
        )

    // Emit mapped / partial / unmapped counts from stats JSON to the log
    PANGENOME_MAP.out.stats
        .map { meta, stats_json ->
            def stats = new groovy.json.JsonSlurper().parse(stats_json.toFile())
            log.info "  Pangenome projection [${meta.id}]: " +
                     "mapped=${stats.mapped ?: '?'}, " +
                     "partial=${stats.partial ?: '?'}, " +
                     "unmapped=${stats.unmapped ?: '?'}"
        }
}
