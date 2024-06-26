import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

task collect_chunks{
    Array[File] input_files
    String output_dir_name

    # Runtime environment
    String docker = "ubuntu:22.04@sha256:a6d2b38300ce017add71440577d5b0a90460d0e57fd7aec21dd0d1b0761bbfb2"
    String ecr = "public.ecr.aws/ubuntu/ubuntu:22.04_stable"
    String container_source = "docker"
    String container_image = if(container_source == "docker") then docker else ecr
    Int cpu = 2
    Int mem_gb = 4
    Int max_retries = 3

    command <<<
        set -e

        # Loop through files in input file and copy/decompress them to output dir
        for input_file in ${sep=" " input_files}; do
            tar -xvzf "$input_file" -C ./
        done

        # Make a list of files in directory
        find ${output_dir_name}/* -type f > ${output_dir_name}.contents.txt

        # Compress directory so it can be placed inside higher-level MultiQC directories
        tar -cvzf ${output_dir_name}.tar.gz ${output_dir_name}
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File output_dir_file_list = "${output_dir_name}.contents.txt"
        File output_dir = "${output_dir_name}.tar.gz"
    }
}

workflow collect_large_file_list_wf{
    Array[File] input_files
    Int chunk_size = 50
    String output_dir_name
    Int num_chunks = ceil(length(input_files)/(chunk_size*1.0))
    Int collect_chunk_cpus = 2
    Int collect_chunk_mem_gb = 4
    String container_source = "docker"

    # Create list of file
    scatter(chunk in range(num_chunks)){

        # Get a slice of the list
        Int start_pos = chunk*chunk_size
        call UTILS.slice{
            input:
                inputs = input_files,
                start_pos = start_pos,
                end_pos = start_pos + chunk_size,
                container_source = container_source
        }

        # zip the sliced files into a tarball
        call UTILS.collect_files as collect_one_chunk{
            input:
                output_dir_name = output_dir_name,
                input_files = slice.outputs,
                container_source = container_source
        }
    }

    # Collect all the scattered chunks into single tarball
    call collect_chunks{
        input:
            input_files = collect_one_chunk.output_dir,
            output_dir_name = output_dir_name,
            cpu = collect_chunk_cpus,
            mem_gb = collect_chunk_mem_gb,
            container_source = container_source
    }

    output{
        File output_dir = collect_chunks.output_dir
    }
}