task make_gwas_summary_stats {
    File file_in_summary_stats
    File file_in_info
    String file_in_summary_stats_format
    String file_in_info_format = "info"
    String file_out_prefix
    Int chunk_size = 50000

    # Optional pop maf files to add pop MAFs from
    File? file_in_pop_mafs
    String? population

    # Runtime options
    String docker = "rtibiocloud/make_gwas_summary_stats:v2_71ab83a"
    Int cpu = 1
    Int mem_gb = 3
    Int max_retries = 3

    command{
        python /opt/make_gwas_summary_stats.py \
            --file_in_summary_stats ${file_in_summary_stats} \
            --file_in_info ${file_in_info} \
            --file_in_summary_stats_format ${file_in_summary_stats_format} \
            --file_in_info_format ${file_in_info_format} \
            ${"--file_in_pop_mafs " + file_in_pop_mafs} \
            ${"--population " + population} \
            --file_out_prefix ${file_out_prefix} \
            --chunk_size ${chunk_size}
    }

    runtime{
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output {
        File output_file = "${file_out_prefix}.tsv.gz"
        File log_file = "${file_out_prefix}.log"
   }


}