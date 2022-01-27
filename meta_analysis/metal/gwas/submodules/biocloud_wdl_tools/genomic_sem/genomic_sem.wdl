task genomic_sem {
    Array[File]? input_files
    Array[String]? trait_names
    Array[String]? sample_sizes
    Array[String]? sample_prev
    Array[String]? pop_prev
    File? reference
    Float? info_filter
    Float? maf_filter
    String out_prefix
    File? ld
    Boolean? munge
    Boolean? LDSC
    File LDSC_file
    String estimation
    Boolean? common_factor
    File? common_factor_model
    Array[String]? se_logit
    Boolean? sumstats
    File sumstats_file
    Boolean common_factor_gwas
    File common_factor_gwas_model
    #Boolean parallel

    # Runtime attributes
    String docker = "rtibiocloud/genomic_sem:v1_75a57cd"
    Int cpu = 16
    Int mem_gb = 10
    Int max_retries = 3

    command {
        set -e
        Rscript /opt/GenomicSEM_commonFactor.R \
            --input_files ${sep="," input_files} \
            --trait_names ${sep="," trait_names} \
            --sample_sizes ${sep="," sample_sizes} \
            --sample_prev ${sep="," sample_prev} \
            --pop_prev ${sep="," pop_prev} \
            --reference ${reference} \
            --info_filter ${info_filter} \
            --maf_filter ${maf_filter} \
            --out_prefix ${out_prefix} \
            --ld ${ld} \
            --munge FALSE \
            --LDSC FALSE \
            --LDSC_file ${LDSC_file} \
            --estimation ${estimation} \
            --common_factor FALSE \
            --common_factor_model ${common_factor_model} \
            --se_logit ${sep="," se_logit} \
            --sumstats FALSE \
            --sumstats_file ${sumstats_file} \
            --common_factor_gwas TRUE \
            --common_factor_gwas_model ${common_factor_gwas_model}
    }

    output {
        File gsem_out = "${out_prefix}.rds"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

}
