task merge{
    # Utility for merging multiple VCF/BCF files
    File file_in

    String output_filename
    String output_type

    # Runtime environment
    String docker = "rtibiocloud/bcftools:v1.9-8875c1e"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        # Merge into single VCF
        bcftools merge \
            --threads ${cpu} \
            -o ${output_filename} \
            -O ${output_type} \
            -l ${file_in}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File vcf_out = "${output_filename}"
    }
}
