import "biocloud_gwas_workflows/association_testing/genomic_sem/genomicSEM_chr_wf.wdl" as GSEM_CHR

workflow genomic_sem_gwas_wf{
    File LDSC_file
    Array[String] out_prefixes
    Array[File] sumstats_files
    Array[String] chrs
    String estimation
    Boolean common_factor_gwas
    File common_factor_gwas_model
    Boolean parallel

    # options
    Array[File]? input_files
    Array[String]? trait_names
    Array[String]? sample_sizes
    Array[String]? sample_prev
    Array[String]? pop_prev
    File? reference
    Float? info_filter
    Float? maf_filter
    File? ld
    Boolean? munge
    Boolean? LDSC
    Boolean? common_factor
    File? common_factor_model
    Array[String]? se_logit

    Int gsem_cpu = 16
    Int gsem_mem_gb = 10

    # Do genomic_sem chr workflow on each chromosome in parallel
    scatter(index in range(length(sumstats_files))){

        call GSEM_CHR.genomicSEM_chr_wf as genomic_sem{
            input:
                input_files = input_files,
                trait_names = trait_names,
                sample_sizes = sample_sizes,
                sample_prev = sample_prev,
                pop_prev = pop_prev,
                reference = reference,
                info_filter = info_filter,
                maf_filter = maf_filter,
                out_prefix = out_prefixes[index],
                ld = ld,
                LDSC_file = LDSC_file,
                estimation = estimation,
                common_factor_model = common_factor_model,
                se_logit = se_logit,
                sumstats_file = sumstats_files[index],
                common_factor_gwas_model = common_factor_gwas_model,
                common_factor_gwas = common_factor_gwas,
                chr = chrs[index],
                parallel = parallel
        }
    }

    output{
        Array[File] genomic_sem_sumstats = genomic_sem.summary_stats
    }
}
