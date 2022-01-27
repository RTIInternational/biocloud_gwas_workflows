task convert_bgen_to_seq_gds {
    File in_bgen
    String out_gds
    String storage_option = "LZMA_RA"
    String float_type = "double"
    Boolean geno = false
    Boolean dosage = false
    Boolean prob = false
    Boolean optimize = false
    Int parallel = 8

    # Runtime environment
    String docker = "rtibiocloud/convert_bgen_to_seq_gds:v1_c3bbaef"
    Int mem_gb = 8
    Int max_retries = 3

    command{
        /opt/convert_bgen_to_seq_gds.R \
            --in-bgen ${in_bgen} \
            --out-gds ${out_gds} \
            --storage-option "${storage_option}" \
            --float-type "${float_type}" \
            ${true="--geno" false="" geno} \
            ${true="--dosage" false="" dosage} \
            ${true="--prob" false="" prob} \
            ${true="--optimize" false="" optimize} \
            --parallel ${parallel}
    }

    runtime{
        docker: docker
        cpu: parallel
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output {
        File output_file = "${out_gds}"
    }
}
