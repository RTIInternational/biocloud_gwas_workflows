task append {
    Array[String] a
    String? b

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        cat ${write_lines(a)}
        ${'echo ' + b}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}

task collect_files{
    Array[File] input_files
    String output_dir_name

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    meta {
        description: "Gather multiple files into a single gzipped tarball that can be unzipped and directly input to MultiQC"
    }

    parameter_meta {
        input_files: "Files to zip"
        output_dir_name: "Name of directory that will be created and tarballed"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_gb: "(optional) the amount of memory (GB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
    }

    command <<<
        set -e

        # Create directory and copy files into directory
        mkdir -p ${output_dir_name}

        # Loop through files in input file and copy/decompress them to output dir
        for input_file in ${sep=" " input_files}; do

            if [[ $input_file == *.tar.gz ]]; then
                # Untar directory into output directory
                tar -xvzf "$input_file" -C ${output_dir_name}
            else
                # Just copy flat files to MultiQC dir
                cp "$input_file" ${output_dir_name}
            fi

        done

        # Make a list of files in directory
        find ${output_dir_name}/* -type f > ${output_dir_name}.contents.txt

        # Compress directory so it can be placed inside higher-level directories
        tar -cvzf ${output_dir_name}.tar.gz ${output_dir_name}

    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File output_dir_file_list = "${output_dir_name}.contents.txt"
        File output_dir = "${output_dir_name}.tar.gz"
    }
}

task slice{
    # Get 0-based slice of an array
    # end_index is exclusive so it works more or less like python
    Array[String] inputs
    Int start_pos
    Int end_pos
    Int slice_size = end_pos - start_pos

    # Make start pos 1-based because of how tail -N + works
    Int actual_start_pos = start_pos + 1

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        tail -n +${actual_start_pos} ${write_lines(inputs)} | head -${slice_size}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[String] outputs = read_lines(stdout())
    }
}

task flatten_string_array {

    Array[Array[String]] array

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
    for line in $(echo ${sep=', ' array}) ; \
    do echo $line | tr -d '"[],' ; done
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[String]  flat_array = read_lines(stdout())
    }
}

task remove_empty_files{
    Array[File] input_files

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        mkdir non_empty_files
        for file in ${sep=' ' input_files}; do
            if [ -s $file ];then
                cp $file non_empty_files
            fi
        done
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[File] non_empty_files = glob("non_empty_files/*")
    }
}

task wc{
    File input_file

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        if [[ ${input_file} =~ \.gz$ ]]
        then
            gunzip -c ${input_file} | wc -l | cut -d" " -f1
        else
            wc -l ${input_file} | cut -d" " -f1
        fi
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Int num_lines = read_int(stdout())
    }
}

task cut{
    File input_file
    String args
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        cut ${args} ${input_file} > ${output_filename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File output_file = "${output_filename}"
    }
}

task get_file_union{
    # Takes a list of files and outputs the union of lines in each file with no duplicates
    Array[File] input_files
    String output_filename

    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        set -e
        for file in ${sep=" " input_files}
        do
            cat $file >> all_files.txt
        done

        # Dedup merged file
        sort all_files.txt | uniq > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File output_file = "${output_filename}"
    }
}

task replace_chr{
    # Replace all occurances of a character in a file with another
    File input_file
    String output_filename
    String char
    String new_char

    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        sed 's/${char}/${new_char}/g' ${input_file} > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File output_file = "${output_filename}"
    }
}

task raise_error{
    # General module for stopping a pipeline with a custom error message
    String msg
    command <<<
        echo '${msg}' > err_msg.txt
        exit 1
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "500 MB"
        maxRetries: 1
    }

    output{
        File err_msg = "err_msg.txt"
    }
}

task cat{
    # Replace all occurances of a character in a file with another
    Array[File] input_files
    String output_filename
    Boolean input_gzipped = false

    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        if [[ '${input_gzipped}' == 'true' ]]; then
            zcat ${sep=" " input_files} > ${output_filename}
        else
            cat ${sep=" " input_files} > ${output_filename}
        fi
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File output_file = "${output_filename}"
    }
}

task get_file_extension{
    # Get a file extension
    # Optionally get more than just the last file extension using fields to determine how many extensions to grab
    String input_file
    Int fields = 1

    command <<<
        echo ${input_file} | awk -F '\.' \
            '{
                extension = ""
                i = NF
                while((i > NF-${fields}) && (substr($i, 1, 1) ~ /[A-Za-z]/)){
                    extension="." $i extension
                    i--
                }
                print extension
            }'
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "100 MB"
    }

    output{
        String extension = read_string(stdout())
    }
}

task array_contains{
    # Return true if array contains an exact match, false otherwise
    Array[String] input_array
    String query

    command<<<
        set -e
        contains=0
        for i in ${sep=" " input_array};do
            if [[ "${query}" == "$i" ]];then
                echo "true"
                contains=1
            fi
        done
        if [ $contains -eq 0 ];then
            echo "false"
        fi
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "100 MB"
    }

    output{
        Boolean contains = read_boolean(stdout())
    }
}

task append_column{
    # Add a additional column to text file where every entry will be 'value'
    File input_file
    String value
    String output_filename

    # Optionally handle different input/output separators
    # Default is whitespace that awk will detect and output using the same separator
    String? OFS
    String? F
    String f_arg = if(defined(F)) then "-F '${F}'" else ""
    String ofs_arg = if(defined(OFS)) then "-v OFS='${OFS}'" else ""

    command<<<
        awk ${f_arg} ${ofs_arg} '{ $(NF+1) = "${value}"; print }' ${input_file} > ${output_filename}
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "100 MB"
    }

    output{
        File output_file = "${output_filename}"
    }
}

task paste{
    Array[File] input_files
    String? delim
    Boolean? s
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        paste ${'-d ' + delim} \
            ${true='-s' false='' s} \
            ${sep=" " input_files} > ${output_filename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File output_file = "${output_filename}"
    }
}

task array_equals{
    # Return true if two arrays contain same values in same order
    Array[String] array_a
    Array[String] array_b

    command<<<
        diff ${write_lines(array_a)} ${write_lines(array_b)} > compare.txt
        if [ -s compare.txt ];then
            echo "false"
        else
            echo "true"
        fi
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "100 MB"
    }

    output{
        Boolean is_equal = read_boolean(stdout())
    }
}

task shuf{
    # Write a random permutation of the input lines to standard output
    File input_file
    Int? n
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command<<<
        shuf ${"-n " + n} ${input_file}  > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output{
        File output_file = "${output_filename}"
    }
}
