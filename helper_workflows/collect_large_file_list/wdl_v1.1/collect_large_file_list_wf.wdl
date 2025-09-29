version 1.1

import "utils.wdl" as UTILS

workflow collect_large_file_list_wf{

    input{
        Array[File] input_files
        Int chunk_size = 50
        String output_dir_name
        Int num_chunks = ceil(length(input_files)/(chunk_size*1.0))
        Int collect_chunk_cpus = 2
        Int collect_chunk_mem_gb = 4
        String image_source = "docker"
        String? ecr_repo
    }

    # Create list of file
    scatter(chunk in range(num_chunks)){

        # Get a slice of the list
        Int start_pos = chunk*chunk_size
        call UTILS.slice{
            input:
                inputs = input_files,
                start_pos = start_pos,
                end_pos = start_pos + chunk_size,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # zip the sliced files into a tarball
        call UTILS.collect_files as collect_one_chunk{
            input:
                output_dir_name = output_dir_name,
                input_files = slice.outputs,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Collect all the scattered chunks into single tarball
    call collect_chunks{
        input:
            input_files = collect_one_chunk.output_dir,
            output_dir_name = output_dir_name,
            cpu = collect_chunk_cpus,
            mem_gb = collect_chunk_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File output_dir = collect_chunks.output_dir
    }
}

task collect_chunks{

    input{
        Array[File] input_files
        String output_dir_name

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:a6d2b38300ce017add71440577d5b0a90460d0e57fd7aec21dd0d1b0761bbfb2"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String image_source = "docker"
        String? ecr_repo
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 2
        Int mem_gb = 4
    }

    command <<<
        set -e

        # Loop through files in input file and copy/decompress them to output dir
        for input_file in ~{sep(" ", input_files)}; do
            tar -xvzf "$input_file" -C ./
        done

        # Make a list of files in directory
        find ~{output_dir_name}/* -type f > ~{output_dir_name}.contents.txt

        # Compress directory so it can be placed inside higher-level MultiQC directories
        tar -cvzf ~{output_dir_name}.tar.gz ~{output_dir_name}
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output{
        File output_dir_file_list = "~{output_dir_name}.contents.txt"
        File output_dir = "~{output_dir_name}.tar.gz"
    }
}
