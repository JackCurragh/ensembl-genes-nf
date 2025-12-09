nextflow.enable.dsl=2
import java.nio.file.Files

include {CHECK_REPEAT_LIBRARY} from '../modules/check_repeat_library.nf'
include {DOWNLOAD_REPEAT_LIBRARY} from '../modules/download_repeat_library.nf'
include {FETCH_GENOME} from '../modules/fetch_genome.nf'
include {GENERATE_REPEATMODELER_LIBRARY} from '../modules/generate_repeatmodeler_library.nf'

workflow {
    // Check if outDir parameter is defined
    if (!params.outDir) {
        error "Undefined --outDir parameter. Please provide the output directory's path"
    }

    // Check if csvFile parameter is defined and exists
    if (params.csvFile) {
        csvFile = file(params.csvFile, checkIfExists: true)
    } else {
        error 'csvFile file not specified!'
    }

    // Read data from the CSV file, split it, and map each row to extract species_name and GCA values
    data_with_gca = Channel.fromPath(csvFile, type: 'file')
           .splitCsv(sep: ',', header: false)
           .map { row -> tuple(row[0], row[1]) }

    // Check if repeat library exists on Ensembl FTP
    library_check = CHECK_REPEAT_LIBRARY(data_with_gca)

    // Branch: separate genomes with existing libraries vs those that need generation
    library_check.branch {
        exists: it[2].trim() == 'exists'
        missing: it[2].trim() == 'missing'
    }.set { branched }

    // Path 1: Download pre-existing libraries (fast, no genome needed)
    existing_libraries = branched.exists
        .map { species, gca, status -> tuple(species, gca) }
    DOWNLOAD_REPEAT_LIBRARY(existing_libraries)

    // Path 2: Fetch genome and generate RepeatModeler library (slow)
    genomes_to_process = branched.missing
        .map { species, gca, status -> tuple(species, gca) }

    fetched_genome = FETCH_GENOME(genomes_to_process)

    GENERATE_REPEATMODELER_LIBRARY(fetched_genome)
}    
