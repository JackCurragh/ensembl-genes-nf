#!/usr/bin/env nextflow
/*
========================================================================================
    VALIDATE_MODELS PIPELINE
========================================================================================
    Score gene models from multiple annotation sources before consolidation.

    For each input GFF3 the pipeline applies:
      1. Structural sanity checks:
         - Splice-site canonicity (GT-AG / GC-AG / AT-AC vs non-canonical)
         - CDS integrity (no in-frame stops, length divisible by 3)
         - ORF completeness (start ATG + stop codon present)
      2. Splice-junction support (when --sj_tabs provided):
         - Fraction of predicted introns backed by ≥1 STAR uniquely-mapped read
         - Median unique-read depth across supported junctions
      3. Optional DIAMOND blastp vs UniProt reviewed (when --protein_db provided):
         - Boosts validation_score by ≤0.1 for transcripts with a strong protein hit

    Composite score (0–1) stored as GFF3 attribute:
      validation_score = 0.4 * structural_score + 0.6 * splice_support_score
      (if no RNA-seq SJ tabs: validation_score = structural_score)

    Inputs:
      --gff3_files     Comma-separated absolute paths to annotation GFF3 files
      --genome_fasta   Softmasked genome FASTA (must have a .fai index)
      --sj_tabs        (optional) Comma-sep STAR SJ.out.tab paths for splice support
      --protein_db     (optional) UniProt reviewed FASTA for DIAMOND blastp
      --outdir         Output directory

    Run locally (stub):
      nextflow run pipelines/validate_models -profile local -stub \
        --gff3_files   'a.gff3,b.gff3' \
        --genome_fasta genome.fa \
        --outdir       results/validate_models
========================================================================================
*/

nextflow.enable.dsl = 2

include { VALIDATE_MODELS } from './modules/validate_models.nf'
include { DIAMOND_BLASTP  } from './modules/diamond_blastp.nf'
include { WRITE_MANIFEST  } from './modules/write_manifest.nf'

params.gff3_files   = null
params.genome_fasta = null
params.sj_tabs      = null   // optional
params.protein_db   = null   // optional; skip DIAMOND if absent

def validate_params() {
    def errors = []
    if (!params.gff3_files)   errors << "  --gff3_files is required"
    if (!params.genome_fasta) errors << "  --genome_fasta is required"
    if (!params.outdir)       errors << "  --outdir is required"
    if (errors) {
        log.error "Missing required parameters:\n${errors.join('\n')}"
        System.exit(1)
    }
}

workflow {

    validate_params()

    // Build (meta, gff3) channel from comma-separated file list.
    // Label = filename without path and without trailing .scored.gff3 / .gff3
    ch_gff3 = Channel
        .fromList(
            params.gff3_files.split(',').collect { it.trim() }.findAll { it }
        )
        .map { p ->
            def label = p.tokenize('/')[-1].replaceAll(/\.scored\.gff3$|\.gff3$/, '')
            tuple([id: label], file(p, checkIfExists: true))
        }

    ch_genome = file(params.genome_fasta, checkIfExists: true)

    // Optional: STAR splice-junction tabs
    def sj_list = params.sj_tabs
        ? params.sj_tabs.split(',').collect { file(it.trim()) }
        : []

    // ── Step 1: structural + splice scoring ──────────────────────────────────
    VALIDATE_MODELS(ch_gff3, ch_genome, sj_list)

    // ── Step 2: DIAMOND protein validation (optional) ────────────────────────
    if (params.protein_db) {
        DIAMOND_BLASTP(
            VALIDATE_MODELS.out.gff3,
            ch_genome,
            file(params.protein_db, checkIfExists: true)
        )
        ch_final = DIAMOND_BLASTP.out.gff3
    } else {
        ch_final = VALIDATE_MODELS.out.gff3
    }

    // ── Step 3: Publish scored GFF3s and write manifest ──────────────────────
    // Collect all (meta, gff3) pairs → pass gff3 paths to WRITE_MANIFEST.
    WRITE_MANIFEST(
        params.outdir,
        ch_final.map { meta, gff3 -> gff3 }.collect()
    )

    // Software versions
    Channel.empty()
        .mix(VALIDATE_MODELS.out.versions)
        .collectFile(
            name:     'software_versions.tsv',
            newLine:  true,
            storeDir: "${params.outdir}/pipeline_info"
        )
}
