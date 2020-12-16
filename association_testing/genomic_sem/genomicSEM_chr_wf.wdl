import "biocloud_gwas_workflows/biocloud_wdl_tools/genomic_sem/genomic_sem.wdl" as GSEM

workflow genomicSEM_chr_wf{
    File LDSC_file
    String out_prefix
    File sumstats_file
    String estimation
    Boolean common_factor_gwas
    File common_factor_gwas_model

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
    String chr              # Chr being analyzed
    Int gsem_cpu = 16
    Int gsem_mem_gb = 10


    call GSEM.genomic_sem{
        input:
            input_files = input_files,
            trait_names = trait_names,
            sample_sizes = sample_sizes,
            sample_prev = sample_prev,
            pop_prev = pop_prev,
            reference = reference,
            info_file = info_filter,
            maf_filter = maf_filter,
            out_prefix = out_prefix,
            ld = ld,
            LDSC_file = LDSC_file,
            estimation = estimation,
            common_factor_model = common_factor_model,
            se_logit = se_logit,
            sumstats_file = sumstats_file,
            common_factor_gwas_model = common_factor_gwas_model,
            common_factor_gwas = common_factor_gwas
    }
    

    output {
        File summary_stats = genomic_sem.sumstats_out
    }
}
