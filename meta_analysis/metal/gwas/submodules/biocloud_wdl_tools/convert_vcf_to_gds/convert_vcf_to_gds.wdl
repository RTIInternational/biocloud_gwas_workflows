task convert_vcf_to_gds {
    File in_file
    String out_file

    # Runtime environment
    String docker = "rtibiocloud/convert_vcf_to_gds:v1_276fc5b"
    Int cpu = 1
    Int mem_gb = 8
    Int max_retries = 3

    command{
        /opt/convert_vcf_to_gds.R \
            --in ${in_file} \
            --out ${out_file}
    }

    runtime{
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output {
        File output_file = "${out_file}"
    }
}
