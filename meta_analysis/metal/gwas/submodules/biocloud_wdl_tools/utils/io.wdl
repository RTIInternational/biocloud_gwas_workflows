task gzip{
    File input_file
    String? user_filename
    String default_filename = basename(input_file) + ".gz"
    String output_filename = select_first([user_filename, default_filename])

    String docker = "rtibiocloud/pigz:v2.4_b243f9"
    Int cpu = 1
    Int mem_gb = 1
    Int max_retries = 3

    command <<<
        pigz -ck -p${cpu} ${input_file} > ${output_filename}
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

task gunzip{
    File input_file
    String? user_filename
    String default_filename = basename(input_file, ".gz")
    String output_filename = select_first([user_filename, default_filename])

    String docker = "rtibiocloud/pigz:v2.4_b243f9"
    Int cpu = 1
    Int mem_gb = 1
    Int max_retries = 3

    command <<<
        unpigz -ck -p${cpu} ${input_file} > ${output_filename}
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
