task tsv_append{
    # TSV utility for concatenating multiple TSVs into one TSV while taking header into account
    File tsv_inputs_tarball
    String output_filename
    Boolean header = true
    Boolean track_source = false
    String? source_header
    String? delimiter
    String tsv_dir = basename(tsv_inputs_tarball, ".tar.gz")


    # Runtime environment
    String docker = "rtibiocloud/tsv-utils:v2.2.0_5141a72"
    Int cpu = 2
    Int mem_gb = 4
    Int max_retries = 3

    command <<<
        # Unzip/decompress files to working directory
        tar -xvzf ${tsv_inputs_tarball} -C ./

        # Unzip any/all gzipped files
        find ${tsv_dir}/ -name '*.gz' | while read file
        do
            echo "Unzipping $file"
            gunzip $file
        done

        # Concat all files together
        tsv-append \
            ${"--source-header " + source_header} \
		    ${true="--header" false="" header} \
		    ${true="--track-source" false="" track_source} \
		    ${"--delimiter '" + delimiter + "'"} \
            ${tsv_dir}/* > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File tsv_output = output_filename
    }
}

task tsv_filter{
    # TSV utility for filtering on multiple columns
    File tsv_input
    String output_filename
    Boolean header = true
    Boolean? invert
    String? delimiter
    Boolean? or_filter  # Evaluate tests as an OR rather than an AND clause.

    # Filtering criteria (see tsv-filter options for more details
    # There are too many options here to parameterize so just pass as a string
    String filter_string

    # Runtime environment
    String docker = "rtibiocloud/tsv-utils:v2.2.0_5141a72"
    Int cpu = 2
    Int mem_gb = 4
    Int max_retries = 3

    command <<<

        input_file=${tsv_input}

        # Unzip tsv input file if necessary
        if [[ ${tsv_input} =~ \.gz$ ]]; then
            log_info "${tsv_input} is gzipped. Unzipping..."
            gunzip -c ${tsv_input} > input.txt
            input_file=input.txt
        fi

        tsv-filter \
		    ${true="--header" false="" header} \
		    ${true="--or" false="" or_filter} \
		    ${true="--invert" false="" invert} \
		    ${"--delimiter '" + delimiter + "'"} \
		    ${filter_string} \
            $input_file > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File tsv_output = output_filename
    }
}

task tsv_select{
    # TSV utility for selecting/re-ordering columns (similar to cut but allows re-ordering)
    File tsv_input
    String output_filename
    Array[String] fields
    Boolean header = true
    String? delimiter
    String rest = "none"  # Location for remaining fields (none|first|last)


    # Runtime environment
    String docker = "rtibiocloud/tsv-utils:v2.2.0_5141a72"
    Int cpu = 2
    Int mem_gb = 4
    Int max_retries = 3

    command <<<

        input_file=${tsv_input}

        # Unzip tsv input file if necessary
        if [[ ${tsv_input} =~ \.gz$ ]]; then
            log_info "${tsv_input} is gzipped. Unzipping..."
            gunzip -c ${tsv_input} > input.txt
            input_file=input.txt
        fi
        tsv-select \
		    ${true="--header" false="" header} \
		    ${"--delimiter '" + delimiter + "'"} \
		    --fields ${sep="," fields} \
		    --rest ${rest} \
            $input_file > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File tsv_output = output_filename
    }
}

task tsv_join{
    # TSV utility for subsetting tsv file based on ids in another file
    File tsv_input
    File tsv_filter_file
    String? key_fields
    String? data_fields
    String? append_fields
    Boolean header = true
    String? delimiter
    String prefix = "none"
    Boolean? write_unmatched
    String? write_unmatched_str
    Boolean? exclude
    Boolean? allow_duplicate_keys
    String output_filename

    # Runtime environment
    String docker = "rtibiocloud/tsv-utils:v2.2.0_5141a72"
    Int cpu = 2
    Int mem_gb = ceil(size(tsv_filter_file, "GiB") * 2) + 2
    Int max_retries = 3

    command <<<

        input_file=${tsv_input}
        filter_file=${tsv_filter_file}

        # Unzip tsv input file if necessary
        if [[ ${tsv_input} =~ \.gz$ ]]; then
            log_info "${tsv_input} is gzipped. Unzipping..."
            gunzip -c ${tsv_input} > input.txt
            input_file=input.txt
        fi

        # Unzip filter file if necessary
        if [[ ${tsv_filter_file} =~ \.gz$ ]]; then
            log_info "${tsv_filter_file} is gzipped. Unzipping..."
            gunzip -c ${tsv_filter_file} > filter_input.txt
            filter_file=filter_input.txt
        fi


        tsv-join \
            --filter-file $filter_file \
            ${"--key-fields " + key_fields} \
            ${"--data-fields " + data_fields} \
            ${"--append-fields " + append_fields} \
		    ${true="--header" false="" header} \
		    ${"--prefix " + prefix} \
		    ${"--delimiter '" + delimiter + "'"} \
            ${true="--write-all " false="" write_unmatched}${write_unmatched_str} \
		    ${true="--exclude" false="" exclude} \
		    ${true="--allow-duplicate-keys" false="" allow_duplicate_keys} \
            $input_file > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File tsv_output = output_filename
    }
}

