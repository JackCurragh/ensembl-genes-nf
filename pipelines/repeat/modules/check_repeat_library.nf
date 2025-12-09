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

process CHECK_REPEAT_LIBRARY {
    tag "$gca:check_library"
    label 'fetch_file'

    input:
    tuple val(species_name), val(gca)

    output:
    tuple val(species_name), val(gca), stdout

    script:
    """
    # Normalize species name: capitalize first letter of genus, lowercase species epithet
    SPECIES_NORMALIZED=\$(echo "${species_name}" | awk -F'_' '{print toupper(substr(\$1,1,1)) tolower(substr(\$1,2)) "_" tolower(\$2)}')
    REPEAT_URL="${params.repeats_ftp_base}/\${SPECIES_NORMALIZED}/${gca}.repeatmodeler.fa"

    if curl --head --silent --fail "\$REPEAT_URL" > /dev/null 2>&1; then
        echo -n "exists"
    else
        echo -n "missing"
    fi
    """
}
