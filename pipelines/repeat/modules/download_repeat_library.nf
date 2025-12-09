#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

process DOWNLOAD_REPEAT_LIBRARY {
    tag "$gca:download_library"
    label 'fetch_file'
    publishDir "${params.outDir}/${gca}/library/", mode: 'copy'

    input:
    tuple val(species_name), val(gca)

    output:
    tuple val(gca), val(species_name), path("${gca}.repeatmodeler.fa")

    script:
    """
    # Normalize species name: capitalize first letter of genus, lowercase species epithet
    SPECIES_NORMALIZED=\$(echo "${species_name}" | awk -F'_' '{print toupper(substr(\$1,1,1)) tolower(substr(\$1,2)) "_" tolower(\$2)}')
    REPEAT_URL="${params.repeats_ftp_base}/\${SPECIES_NORMALIZED}/${gca}.repeatmodeler.fa"
    curl --silent --fail --output "${gca}.repeatmodeler.fa" "\$REPEAT_URL"
    echo "Successfully downloaded RepeatModeler library for ${gca}"
    """
}
