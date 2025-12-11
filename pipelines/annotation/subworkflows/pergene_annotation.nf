/*
 * PER-GENE ANNOTATION SUBWORKFLOW
 *
 * Gene-specific translon annotation using profile clustering and custom methods.
 *
 * WORKFLOW:
 * 1. Load count matrix (unique RPFs × samples)
 * 2. For each gene, cluster samples by translation profile similarity
 * 3. Extract gene-specific BAM slices for each cluster
 * 4. Generate cluster-specific translation profiles
 * 5. Run custom profile-based translon detection methods
 * 6. Run statistical tests for differential translation
 * 7. Score and filter translons by cluster support
 * 8. Generate per-gene, per-cluster annotations
 *
 * OUTPUT: Multiple annotations per gene (one per cluster)
 *         Each gene can have different translons in different conditions
 */

include { LOAD_COUNT_MATRIX } from '../modules/pergene/load_count_matrix.nf'
include { CLUSTER_PERGENE } from '../modules/pergene/cluster_pergene.nf'
include { EXTRACT_GENE_BAMS } from '../modules/pergene/extract_gene_bams.nf'
include { GENERATE_CLUSTER_PROFILES } from '../modules/pergene/generate_cluster_profiles.nf'
include { PROFILE_BASED_DETECTION } from '../modules/pergene/profile_based_detection.nf'
include { STATISTICAL_TESTS } from '../modules/pergene/statistical_tests.nf'
include { SCORE_CLUSTER_TRANSLONS } from '../modules/pergene/score_cluster_translons.nf'
include { GENERATE_PERGENE_ANNOTATION } from '../modules/pergene/generate_pergene_annotation.nf'

workflow PERGENE_ANNOTATION {
    take:
    count_matrix            // path: Unique RPFs × samples matrix (parquet)
    merged_bam              // path: All samples merged BAM
    sample_offsets          // path: Per-sample, per-readlength offsets
    reference_gtf           // path: Reference GTF annotation
    reference_fasta         // path: Reference genome FASTA
    rnaseq_results_dir      // path: Optional RNA-Seq results
    similarity_method       // val: 'correlation', 'cosine', 'euclidean'
    similarity_threshold    // val: Per-gene clustering threshold
    min_cluster_size        // val: Minimum samples per cluster

    main:
    // TODO: Step 1 - Load and prepare count matrix
    LOAD_COUNT_MATRIX(
        count_matrix,
        reference_gtf
    )

    // TODO: Step 2 - Cluster samples separately for each gene
    // This produces: gene_id → [cluster_1_samples, cluster_2_samples, ...]
    CLUSTER_PERGENE(
        LOAD_COUNT_MATRIX.out.gene_profiles,
        similarity_method,
        similarity_threshold,
        min_cluster_size
    )

    // TODO: Step 3 - Extract gene-specific BAM slices per cluster
    // For each (gene, cluster) pair, extract reads from samples in that cluster
    EXTRACT_GENE_BAMS(
        merged_bam,
        CLUSTER_PERGENE.out.gene_cluster_assignments,
        reference_gtf
    )

    // TODO: Step 4 - Generate translation profiles for each gene-cluster
    GENERATE_CLUSTER_PROFILES(
        EXTRACT_GENE_BAMS.out.gene_cluster_bams,
        sample_offsets,
        reference_gtf,
        reference_fasta
    )

    // TODO: Step 5 - Profile-based translon detection
    // Custom methods optimized for profile analysis:
    // - Periodicity scoring
    // - Start codon enrichment
    // - Frame preference analysis
    PROFILE_BASED_DETECTION(
        GENERATE_CLUSTER_PROFILES.out.profiles,
        reference_gtf,
        reference_fasta
    )

    // TODO: Step 6 - Statistical tests for translation
    // Compare clusters to find differential translation:
    // - Coverage differences
    // - Frame usage differences
    // - Start codon usage
    STATISTICAL_TESTS(
        GENERATE_CLUSTER_PROFILES.out.profiles,
        CLUSTER_PERGENE.out.gene_cluster_assignments,
        reference_gtf
    )

    // TODO: Step 7 - Score and filter translons
    // Combine evidence from profile analysis and statistical tests
    SCORE_CLUSTER_TRANSLONS(
        PROFILE_BASED_DETECTION.out.candidate_translons,
        STATISTICAL_TESTS.out.significance_scores,
        CLUSTER_PERGENE.out.gene_cluster_assignments
    )

    // TODO: Step 8 - Generate final per-gene annotations
    // Output format: gene_id, cluster_id, translons, evidence
    GENERATE_PERGENE_ANNOTATION(
        SCORE_CLUSTER_TRANSLONS.out.filtered_translons,
        CLUSTER_PERGENE.out.gene_cluster_assignments,
        reference_gtf,
        reference_fasta
    )

    emit:
    gene_cluster_assignments = CLUSTER_PERGENE.out.gene_cluster_assignments
    cluster_profiles = GENERATE_CLUSTER_PROFILES.out.profiles
    pergene_annotations = GENERATE_PERGENE_ANNOTATION.out.annotations
    cluster_stats = CLUSTER_PERGENE.out.stats
    detection_stats = PROFILE_BASED_DETECTION.out.stats
}
