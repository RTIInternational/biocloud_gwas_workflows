task convert_variant_ids {
    File in_file
    File ref
    String chr
    Int in_header = 0
    String in_sep = "tab"
    Int in_id_col
    Int in_chr_col
    Int in_pos_col
    Int in_a1_col
    Int in_a2_col
    String in_missing_allele
    String in_deletion_allele
    String ref_deletion_allele
    Int in_chunk_size = 50000
    Int ref_chunk_size = 1000000
    Boolean? rescue_rsids

    String? output_compression
    String output_filename
    String log_filename = "${output_filename}.log"

    # Runtime environment
    String docker = "rtibiocloud/convert_variant_ids:v1_9a23978"
    Int cpu = 1
    Int mem_gb = 1
    Int max_retries = 3

    command{
        python /opt/convert_variant_ids.py \
            --chr ${chr} \
            --in_file ${in_file} \
            --in_header ${in_header} \
            --in_sep ${in_sep} \
            --in_id_col ${in_id_col} \
            --in_chr_col ${in_chr_col} \
            --in_pos_col ${in_pos_col} \
            --in_a1_col ${in_a1_col} \
            --in_a2_col ${in_a2_col} \
            --in_missing_allele "${in_missing_allele}" \
            --in_deletion_allele "${in_deletion_allele}" \
            --in_chunk_size ${in_chunk_size} \
            --ref ${ref} \
            --ref_deletion_allele "${ref_deletion_allele}" \
            --ref_chunk_size ${ref_chunk_size} \
            ${'--out_compression ' + output_compression} \
            ${true="--rescue_rsids" false="" rescue_rsids} \
            --out_file ${output_filename} \
            --log_file ${log_filename}
    }

    runtime{
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output {
        File output_file = "${output_filename}"
        File log = "${log_filename}"
    }
}
