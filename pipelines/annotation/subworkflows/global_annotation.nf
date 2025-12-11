/*
 * GLOBAL ANNOTATION SUBWORKFLOW
 *
 * Standard translon annotation using published tools and global sample clustering.
 *
 * WORKFLOW:
 * 1. Load Ribo-Seq results (BAMs, offsets, metadata)
 * 2. Compute global translation profiles across all samples
 * 3. Cluster samples by similarity → N meta-samples
 * 4. Merge BAMs and offsets within each cluster
 * 5. Optional: Filter transcripts using RNA-Seq data
 * 6. Run published translon callers on each cluster
 * 7. Standardize predictions to common format
 * 8. Generate consensus annotation per cluster
 *
 * OUTPUT: Standard genome-wide translon annotations
 */

include { LOAD_RIBOSEQ_DATA } from '../modules/global/load_riboseq_data.nf'
include { COMPUTE_GLOBAL_PROFILES } from '../modules/global/compute_global_profiles.nf'
include { CLUSTER_SAMPLES_GLOBAL } from '../modules/global/cluster_samples_global.nf'
include { MERGE_CLUSTER_BAMS } from '../modules/global/merge_cluster_bams.nf'
include { AGGREGATE_CLUSTER_OFFSETS } from '../modules/global/aggregate_cluster_offsets.nf'
include { FILTER_TRANSCRIPTS } from '../modules/global/filter_transcripts.nf'
include { CALL_TRANSLONS } from '../modules/global/call_translons.nf'
include { CONSENSUS_ANNOTATION } from '../modules/global/consensus_annotation.nf'

workflow GLOBAL_ANNOTATION {
    take:
    riboseq_results_dir     // path: Ribo-Seq pipeline results
    reference_gtf           // path: Reference GTF annotation
    reference_fasta         // path: Reference genome FASTA
    rnaseq_results_dir      // path: Optional RNA-Seq results
    similarity_method       // val: 'correlation', 'cosine', 'euclidean'
    similarity_threshold    // val: Clustering threshold (0.0-1.0)
    translon_callers        // val: Comma-separated list of tools
    min_caller_agreement    // val: Minimum callers for consensus

    main:
    // TODO: Step 1 - Load all Ribo-Seq samples
    LOAD_RIBOSEQ_DATA(riboseq_results_dir)

    // TODO: Step 2 - Compute translation profiles for clustering
    COMPUTE_GLOBAL_PROFILES(
        LOAD_RIBOSEQ_DATA.out.bams,
        LOAD_RIBOSEQ_DATA.out.offsets,
        reference_gtf
    )

    // TODO: Step 3 - Cluster samples globally
    CLUSTER_SAMPLES_GLOBAL(
        COMPUTE_GLOBAL_PROFILES.out.profiles,
        LOAD_RIBOSEQ_DATA.out.metadata,
        similarity_method,
        similarity_threshold
    )

    // TODO: Step 4 - Merge BAMs within each cluster
    MERGE_CLUSTER_BAMS(
        LOAD_RIBOSEQ_DATA.out.bams,
        CLUSTER_SAMPLES_GLOBAL.out.cluster_assignments
    )

    // TODO: Step 5 - Aggregate offsets per cluster
    AGGREGATE_CLUSTER_OFFSETS(
        LOAD_RIBOSEQ_DATA.out.offsets,
        CLUSTER_SAMPLES_GLOBAL.out.cluster_assignments,
        LOAD_RIBOSEQ_DATA.out.read_depths
    )

    // TODO: Step 6 - Filter transcripts (optional RNA-Seq enhancement)
    FILTER_TRANSCRIPTS(
        reference_gtf,
        rnaseq_results_dir
    )

    // TODO: Step 7 - Call translons with published tools
    CALL_TRANSLONS(
        MERGE_CLUSTER_BAMS.out.merged_bams,
        AGGREGATE_CLUSTER_OFFSETS.out.cluster_offsets,
        FILTER_TRANSCRIPTS.out.filtered_gtf,
        reference_fasta,
        rnaseq_results_dir,
        translon_callers
    )

    // TODO: Step 8 - Generate consensus annotation
    CONSENSUS_ANNOTATION(
        CALL_TRANSLONS.out.standardized_predictions,
        FILTER_TRANSCRIPTS.out.filtered_gtf,
        reference_fasta,
        min_caller_agreement
    )

    emit:
    cluster_assignments = CLUSTER_SAMPLES_GLOBAL.out.cluster_assignments
    merged_bams = MERGE_CLUSTER_BAMS.out.merged_bams
    final_annotations = CONSENSUS_ANNOTATION.out.final_gtf
    annotation_stats = CONSENSUS_ANNOTATION.out.stats
}
