task convert_seq_gds_to_snp_gds {
    File in_seq_gds
    String out_snp_gds
    String compress_geno = "ZIP_RA"
    Boolean optimize = false
    Boolean dosage = false

    # Runtime environment
    String docker = "rtibiocloud/convert_seq_gds_to_snp_gds:v1_3e0e0c9"
    Int mem_gb = 8
    Int cpu = 1
    Int max_retries = 3

    command{
        /opt/convert_seq_gds_to_snp_gds.R \
            --in-seq-gds ${in_seq_gds} \
            --out-snp-gds ${out_snp_gds} \
            --compress-geno "${compress_geno}" \
            ${true="--optimize" false="" optimize} \
            ${true="--dosage" false="" dosage}
    }

    runtime{
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output {
        File output_file = "${out_snp_gds}"
    }
}
