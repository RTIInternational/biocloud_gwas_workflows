task split_by_variant{
    # Utility for providing variant lists for splitting a GDS file into chunks of N variants.
    File input_gds
    Int chunk_size
    String variant_id_field
    String output_basename

    # Runtime environment
    String docker = "rtibiocloud/split_gds_by_variant:v3.11_ee7a5e8"
    Int cpu = 1
    Int mem_gb = 1
    Int max_retries = 3

    command {
        /opt/split_gds_by_variant.R \
            --file-gds ${input_gds} \
            --variant-id-field ${variant_id_field} \
            --chunk-size ${chunk_size} \
            --out-prefix ${output_basename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        Array[File] split_lists = glob("${output_basename}*")
    }
}
