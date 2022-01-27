task bgzip{
    File input_file
    String? user_filename
    String default_filename = basename(input_file) + ".gz"
    String output_filename = select_first([user_filename, default_filename])

    String docker = "rtibiocloud/htslib:v1.9_cc494f9"
    Int cpu = 1
    Int mem_gb = 1
    Int max_retries = 3

    command <<<
        bgzip -@ ${cpu} -c ${input_file} > ${output_filename}
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

task bgunzip{
    File input_file
    String? user_filename
    String default_filename = basename(input_file, ".gz")
    String output_filename = select_first([user_filename, default_filename])

    String docker = "rtibiocloud/htslib:v1.9_cc494f9"
    Int cpu = 1
    Int mem_gb = 1
    Int max_retries = 3

    command <<<
        bgzip -@ ${cpu} -d -c ${input_file} > ${output_filename}
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