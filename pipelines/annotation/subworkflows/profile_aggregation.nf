/*
 * PROFILE AGGREGATION SUBWORKFLOW
 *
 * Clusters Ribo-Seq samples by translation profile similarity and aggregates
 * similar samples to increase signal-to-noise for ORF detection.
 *
 * Steps:
 * 1. Load Ribo-Seq BAMs and offsets from results directory
 * 2. Compute translation profiles for each sample
 * 3. Calculate pairwise similarity between profiles
 * 4. Cluster samples by similarity threshold
 * 5. Merge BAMs within each cluster
 * 6. Generate aggregated offsets for merged samples
 */

include { LOAD_RIBOSEQ_DATA } from '../modules/load_riboseq_data.nf'
include { COMPUTE_TRANSLATION_PROFILES } from '../modules/compute_translation_profiles.nf'
include { CALCULATE_PROFILE_SIMILARITY } from '../modules/calculate_profile_similarity.nf'
include { CLUSTER_SAMPLES } from '../modules/cluster_samples.nf'
include { MERGE_SAMPLE_BAMS } from '../modules/merge_sample_bams.nf'
include { AGGREGATE_OFFSETS } from '../modules/aggregate_offsets.nf'

workflow PROFILE_AGGREGATION {
    take:
    riboseq_results_dir     // path: Ribo-Seq results directory
    filtered_gtf            // path: Filtered GTF from transcript filtering
    similarity_method       // val: 'correlation', 'cosine', or 'euclidean'
    similarity_threshold    // val: Threshold for clustering (0.0-1.0)

    main:
    // Load Ribo-Seq BAMs and offsets
    LOAD_RIBOSEQ_DATA(riboseq_results_dir)

    riboseq_bams = LOAD_RIBOSEQ_DATA.out.bams
    riboseq_offsets = LOAD_RIBOSEQ_DATA.out.offsets
    sample_metadata = LOAD_RIBOSEQ_DATA.out.metadata

    // Compute translation profiles for each sample
    // Profiles capture codon-level translation patterns across transcriptome
    COMPUTE_TRANSLATION_PROFILES(
        riboseq_bams.combine(riboseq_offsets, by: 0),
        filtered_gtf
    )

    // Calculate pairwise similarity between all samples
    all_profiles = COMPUTE_TRANSLATION_PROFILES.out.profiles
        .map { meta, profile -> profile }
        .collect()

    CALCULATE_PROFILE_SIMILARITY(
        all_profiles,
        similarity_method
    )

    // Cluster samples by similarity threshold
    CLUSTER_SAMPLES(
        CALCULATE_PROFILE_SIMILARITY.out.similarity_matrix,
        sample_metadata,
        similarity_threshold
    )

    // Group BAMs by cluster assignment
    bams_by_cluster = riboseq_bams
        .combine(CLUSTER_SAMPLES.out.sample_assignments)
        .filter { bam_meta, bam, bai, assign_meta, cluster_id ->
            bam_meta.id == assign_meta.sample_id
        }
        .map { bam_meta, bam, bai, assign_meta, cluster_id ->
            [cluster_id, bam_meta, bam, bai]
        }
        .groupTuple()

    // Merge BAMs within each cluster
    MERGE_SAMPLE_BAMS(bams_by_cluster)

    // Aggregate offsets for each cluster (weighted average by read depth)
    offsets_by_cluster = riboseq_offsets
        .combine(CLUSTER_SAMPLES.out.sample_assignments)
        .filter { offset_meta, offsets, assign_meta, cluster_id ->
            offset_meta.id == assign_meta.sample_id
        }
        .map { offset_meta, offsets, assign_meta, cluster_id ->
            [cluster_id, offset_meta, offsets]
        }
        .groupTuple()

    AGGREGATE_OFFSETS(
        offsets_by_cluster,
        LOAD_RIBOSEQ_DATA.out.read_depths
    )

    emit:
    aggregated_bams = MERGE_SAMPLE_BAMS.out.merged_bams          // tuple: [ cluster_id, bam, bai ]
    aggregated_offsets = AGGREGATE_OFFSETS.out.aggregated_offsets // tuple: [ cluster_id, offsets ]
    sample_groupings = CLUSTER_SAMPLES.out.cluster_assignments    // path: cluster_assignments.tsv
    similarity_matrix = CALCULATE_PROFILE_SIMILARITY.out.similarity_matrix
    clustering_stats = CLUSTER_SAMPLES.out.stats
}
