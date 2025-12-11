#!/usr/bin/env nextflow
/*
========================================================================================
    ANNOTATION PIPELINE
========================================================================================
    Translon annotation workflow with two approaches:

    1. GLOBAL ANNOTATION: Sample clustering → published tools → consensus
       - Clusters samples globally into N meta-samples
       - Runs established tools (ORFquant, RiboCode, Ribotish, etc.)
       - Standard genome-wide translon annotation per cluster

    2. PER-GENE ANNOTATION: Gene-level clustering → custom methods → profiles
       - Clusters samples per-gene based on translation profiles
       - Multiple annotations per gene (one per cluster)
       - Custom profile-based and statistical methods
       - Finds condition-specific translons
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

// Validate parameters against schema
validateParameters()

// Check for required parameters
if (!params.gtf) {
    error "Reference GTF annotation file (--gtf) is required"
}
if (!params.fasta) {
    error "Reference genome FASTA file (--fasta) is required"
}

/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

// Global annotation branch
include { GLOBAL_ANNOTATION } from './subworkflows/global_annotation.nf'

// Per-gene annotation branch
include { PERGENE_ANNOTATION } from './subworkflows/pergene_annotation.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // ========================================================================
    // BRANCH 1: GLOBAL ANNOTATION
    // ========================================================================
    // Standard approach compatible with published tools
    // Input: Ribo-Seq results from upstream pipeline
    // Output: Genome-wide translon annotations (one per cluster)

    if (params.mode in ['global', 'both']) {
        GLOBAL_ANNOTATION(
            file(params.riboseq_results_dir),
            file(params.gtf),
            file(params.fasta),
            params.rnaseq_results_dir ? file(params.rnaseq_results_dir) : null,
            params.similarity_method ?: 'correlation',
            params.similarity_threshold ?: 0.8,
            params.translon_callers ?: 'orfquant,ribocode,ribotish',
            params.min_caller_agreement ?: 2
        )
    }

    // ========================================================================
    // BRANCH 2: PER-GENE ANNOTATION
    // ========================================================================
    // Custom approach for gene-specific translation patterns
    // Input: Count matrix + merged BAM from Ribo-Seq pipeline
    // Output: Multiple annotations per gene (condition-specific)

    if (params.mode in ['pergene', 'both']) {
        PERGENE_ANNOTATION(
            file(params.count_matrix),
            file(params.merged_bam),
            file(params.sample_offsets),
            file(params.gtf),
            file(params.fasta),
            params.rnaseq_results_dir ? file(params.rnaseq_results_dir) : null,
            params.pergene_similarity_method ?: 'correlation',
            params.pergene_similarity_threshold ?: 0.8,
            params.min_cluster_size ?: 3
        )
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
