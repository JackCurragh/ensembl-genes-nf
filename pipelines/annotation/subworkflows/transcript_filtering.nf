/*
 * TRANSCRIPT FILTERING SUBWORKFLOW
 *
 * Filters reference GTF to biologically relevant isoforms using:
 * - RNA-Seq expression levels (if available)
 * - RNA-Seq assembled transcripts (if available)
 * - Splice junction support (if available)
 * - Ribo-Seq translation evidence
 *
 * Resolves isoform ambiguity by selecting most supported transcripts
 */

include { LOAD_RNASEQ_DATA } from '../modules/load_rnaseq_data.nf'
include { FILTER_BY_EXPRESSION } from '../modules/filter_by_expression.nf'
include { MERGE_NOVEL_TRANSCRIPTS } from '../modules/merge_novel_transcripts.nf'
include { RESOLVE_ISOFORMS } from '../modules/resolve_isoforms.nf'
include { ANNOTATE_TRANSCRIPT_SUPPORT } from '../modules/annotate_transcript_support.nf'

workflow TRANSCRIPT_FILTERING {
    take:
    reference_gtf           // path: Reference GTF annotation
    reference_fasta         // path: Reference genome FASTA
    rnaseq_results_dir      // path: Optional RNA-Seq results directory
    min_expression_threshold // val: Minimum expression level (TPM/FPKM)

    main:
    // Load RNA-Seq data if available
    if (rnaseq_results_dir) {
        LOAD_RNASEQ_DATA(rnaseq_results_dir)

        rnaseq_expression = LOAD_RNASEQ_DATA.out.expression
        rnaseq_junctions = LOAD_RNASEQ_DATA.out.junctions
        rnaseq_assembly = LOAD_RNASEQ_DATA.out.assembly
    } else {
        rnaseq_expression = Channel.empty()
        rnaseq_junctions = Channel.empty()
        rnaseq_assembly = Channel.empty()
    }

    // Filter transcripts by expression if RNA-Seq data available
    if (rnaseq_results_dir) {
        FILTER_BY_EXPRESSION(
            reference_gtf,
            rnaseq_expression,
            min_expression_threshold
        )
        expressed_gtf = FILTER_BY_EXPRESSION.out.filtered_gtf
    } else {
        // If no RNA-Seq, use full reference GTF
        expressed_gtf = Channel.value(reference_gtf)
    }

    // Merge novel transcripts from RNA-Seq assembly if available
    if (rnaseq_results_dir) {
        MERGE_NOVEL_TRANSCRIPTS(
            expressed_gtf,
            rnaseq_assembly,
            reference_fasta
        )
        merged_gtf = MERGE_NOVEL_TRANSCRIPTS.out.merged_gtf
    } else {
        merged_gtf = expressed_gtf
    }

    // Resolve isoform ambiguity using all available evidence
    RESOLVE_ISOFORMS(
        merged_gtf,
        rnaseq_expression,
        rnaseq_junctions
    )

    // Annotate transcripts with support evidence
    ANNOTATE_TRANSCRIPT_SUPPORT(
        RESOLVE_ISOFORMS.out.resolved_gtf,
        rnaseq_expression,
        rnaseq_junctions
    )

    emit:
    filtered_gtf = ANNOTATE_TRANSCRIPT_SUPPORT.out.annotated_gtf
    filtering_stats = ANNOTATE_TRANSCRIPT_SUPPORT.out.stats
    transcript_support_table = ANNOTATE_TRANSCRIPT_SUPPORT.out.support_table
}
