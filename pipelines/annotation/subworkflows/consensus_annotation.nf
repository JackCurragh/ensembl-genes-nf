/*
 * CONSENSUS ANNOTATION SUBWORKFLOW
 *
 * Reconciles ORF predictions from multiple callers to generate final annotation:
 *
 * 1. Overlap analysis: Find regions predicted by multiple callers
 * 2. Scoring: Rank ORFs by:
 *    - Number of callers in agreement
 *    - Translation metrics (coverage, periodicity, read density)
 *    - Conservation scores (if provided)
 *    - Consistency with reference annotation
 * 3. Conflict resolution: Handle overlapping/conflicting predictions
 * 4. Annotation generation: Create final GTF with confidence scores
 * 5. Quality metrics: Generate QC reports and statistics
 */

include { OVERLAP_ORF_PREDICTIONS } from '../modules/overlap_orf_predictions.nf'
include { SCORE_ORFS } from '../modules/score_orfs.nf'
include { RESOLVE_CONFLICTS } from '../modules/resolve_conflicts.nf'
include { GENERATE_FINAL_GTF } from '../modules/generate_final_gtf.nf'
include { ANNOTATE_ORF_FEATURES } from '../modules/annotate_orf_features.nf'
include { GENERATE_QC_REPORT } from '../modules/generate_qc_report.nf'

workflow CONSENSUS_ANNOTATION {
    take:
    orf_predictions         // tuple: [ meta, standardized_predictions_by_caller ]
    filtered_gtf            // path: Filtered reference GTF
    reference_fasta         // path: Reference genome FASTA
    min_caller_agreement    // val: Minimum number of callers to support ORF
    conservation_scores     // path: Optional conservation scores (PhyloP, PhastCons)

    main:
    // Find overlapping ORF predictions across callers
    OVERLAP_ORF_PREDICTIONS(
        orf_predictions,
        filtered_gtf
    )

    // Score each ORF based on multiple evidence types
    SCORE_ORFS(
        OVERLAP_ORF_PREDICTIONS.out.overlapping_orfs,
        orf_predictions,
        filtered_gtf,
        conservation_scores
    )

    // Resolve conflicts between overlapping ORFs
    // Priority: reference > high-confidence novel > low-confidence novel
    RESOLVE_CONFLICTS(
        SCORE_ORFS.out.scored_orfs,
        filtered_gtf,
        min_caller_agreement
    )

    // Generate final GTF annotation with confidence scores and metadata
    GENERATE_FINAL_GTF(
        RESOLVE_CONFLICTS.out.resolved_orfs,
        filtered_gtf,
        reference_fasta
    )

    // Annotate ORF features (uORFs, dORFs, overlapping ORFs, novel ORFs)
    ANNOTATE_ORF_FEATURES(
        GENERATE_FINAL_GTF.out.final_gtf,
        filtered_gtf
    )

    // Generate comprehensive QC report
    GENERATE_QC_REPORT(
        ANNOTATE_ORF_FEATURES.out.annotated_gtf,
        OVERLAP_ORF_PREDICTIONS.out.overlap_stats,
        SCORE_ORFS.out.scoring_stats,
        RESOLVE_CONFLICTS.out.conflict_stats,
        orf_predictions
    )

    emit:
    final_annotation = ANNOTATE_ORF_FEATURES.out.annotated_gtf
    orf_scores = SCORE_ORFS.out.scored_orfs
    annotation_stats = GENERATE_QC_REPORT.out.stats
    qc_report = GENERATE_QC_REPORT.out.report
    orf_feature_table = ANNOTATE_ORF_FEATURES.out.feature_table
}
