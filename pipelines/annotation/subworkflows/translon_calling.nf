/*
 * ORF CALLING SUBWORKFLOW
 *
 * Runs multiple ORF prediction tools on aggregated Ribo-Seq data:
 * - ORFquant: Translation-based ORF detection
 * - RiboCode: ORF detection using ribosome profiling data
 * - Ribotish: ORF calling with quality scores
 * - ORF-RATER: Machine learning-based ORF classification
 * - DeepRibo: Deep learning ORF prediction (optional)
 *
 * Each caller receives:
 * - Aggregated Ribo-Seq BAMs
 * - P-site offsets
 * - Filtered transcript GTF
 * - Optional RNA-Seq data for background/expression
 */

include { PREPARE_ORF_CALLING_INPUT } from '../modules/prepare_orf_calling_input.nf'
include { RUN_ORFQUANT } from '../modules/run_orfquant.nf'
include { RUN_RIBOCODE } from '../modules/run_ribocode.nf'
include { RUN_RIBOTISH } from '../modules/run_ribotish.nf'
include { RUN_ORFRATER } from '../modules/run_orfrater.nf'
include { RUN_DEEPRIBO } from '../modules/run_deepribo.nf'
include { STANDARDIZE_ORF_PREDICTIONS } from '../modules/standardize_orf_predictions.nf'

workflow ORF_CALLING {
    take:
    aggregated_bams         // tuple: [ cluster_id, bam, bai ]
    aggregated_offsets      // tuple: [ cluster_id, offsets ]
    filtered_gtf            // path: Filtered GTF annotation
    reference_fasta         // path: Reference genome FASTA
    rnaseq_results_dir      // path: Optional RNA-Seq results for background
    orf_callers             // val: Comma-separated list of callers to run

    main:
    // Parse caller list
    callers_list = orf_callers.split(',').collect { it.trim() }

    // Combine BAMs with offsets by cluster_id
    bam_offset_combined = aggregated_bams
        .combine(aggregated_offsets, by: 0)
        .map { cluster_id, bam, bai, offsets ->
            [[id: cluster_id], bam, bai, offsets]
        }

    // Prepare standardized input for all ORF callers
    PREPARE_ORF_CALLING_INPUT(
        bam_offset_combined,
        filtered_gtf,
        reference_fasta
    )

    // Initialize empty channel for predictions
    all_predictions = Channel.empty()

    // Run ORFquant if requested
    if ('orfquant' in callers_list) {
        RUN_ORFQUANT(
            PREPARE_ORF_CALLING_INPUT.out.prepared_data,
            filtered_gtf,
            reference_fasta
        )
        all_predictions = all_predictions.mix(
            RUN_ORFQUANT.out.predictions.map { meta, pred -> ['orfquant', meta, pred] }
        )
    }

    // Run RiboCode if requested
    if ('ribocode' in callers_list) {
        RUN_RIBOCODE(
            PREPARE_ORF_CALLING_INPUT.out.prepared_data,
            filtered_gtf,
            reference_fasta
        )
        all_predictions = all_predictions.mix(
            RUN_RIBOCODE.out.predictions.map { meta, pred -> ['ribocode', meta, pred] }
        )
    }

    // Run Ribotish if requested
    if ('ribotish' in callers_list) {
        RUN_RIBOTISH(
            PREPARE_ORF_CALLING_INPUT.out.prepared_data,
            filtered_gtf,
            reference_fasta
        )
        all_predictions = all_predictions.mix(
            RUN_RIBOTISH.out.predictions.map { meta, pred -> ['ribotish', meta, pred] }
        )
    }

    // Run ORF-RATER if requested
    if ('orfrater' in callers_list) {
        RUN_ORFRATER(
            PREPARE_ORF_CALLING_INPUT.out.prepared_data,
            filtered_gtf,
            reference_fasta,
            rnaseq_results_dir
        )
        all_predictions = all_predictions.mix(
            RUN_ORFRATER.out.predictions.map { meta, pred -> ['orfrater', meta, pred] }
        )
    }

    // Run DeepRibo if requested
    if ('deepribo' in callers_list) {
        RUN_DEEPRIBO(
            PREPARE_ORF_CALLING_INPUT.out.prepared_data,
            filtered_gtf,
            reference_fasta
        )
        all_predictions = all_predictions.mix(
            RUN_DEEPRIBO.out.predictions.map { meta, pred -> ['deepribo', meta, pred] }
        )
    }

    // Standardize all predictions to common format (BED12/GTF with scores)
    STANDARDIZE_ORF_PREDICTIONS(
        all_predictions.groupTuple(by: 1),  // Group by meta (cluster_id)
        filtered_gtf
    )

    emit:
    orf_predictions = STANDARDIZE_ORF_PREDICTIONS.out.standardized_predictions
    individual_caller_outputs = all_predictions
    calling_stats = STANDARDIZE_ORF_PREDICTIONS.out.stats
}
