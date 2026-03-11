/*
 * Simple Sequential Subworkflow Example
 *
 * This shows a subworkflow where processes run sequentially
 * (output of one becomes input to the next)
 */

include { PREP_AA_FASTA } from '../modules/prep_aa_fasta'
include { PREP_DIAMOND_DB } from '../modules/diamond_db_prep'
include { DIAMOND } from '../modules/diamond'

workflow DIAMOND_PROTEIN_VALIDATION {
    take:
    genome_fa           // channel: [ val(meta), path(files) ]
    annotation_gff3     // channel: [ val(meta), path(files) ]
    ref_protein_faa     // channel: [ val(meta), path(files) ]

    main:
    PREP_AA_FASTA( genome_fa, annotation_gff3)
    PREP_DIAMOND_DB( ref_protein_faa )
    DIAMOND( PREP_AA_FASTA.out.aa_fasta, PREP_DIAMOND_DB.out.db )
    PARSE_DIAMOND_OUTPUT( DIAMOND.out.diamond_output )

    emit:
    report = DIAMOND.out.diamond_output  // channel: [ val(meta), path(combined_file) ]
    versions = DIAMOND.out.versions.concat(PREP_AA_FASTA.out.versions)
}
