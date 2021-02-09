import "biocloud_gwas_workflows/helper_workflows/cov_ldsc/cov_ldsc_reference_chr_wf.wdl" as COVLDSC_CHR

workflow cov_ldsc_wf{
    Array[String] plink_format_prefixes
    File = cov_file
    Array[String] output_file_prefixes
    Array[String] chrs
    Int ldsc_cpu
    Int ldsc_mem_gb


    # Do cov-LDSC chr workflow on each chromosome in parallel
    scatter(index in range(length(plink_format_prefixes))){

        call COVLDSC_CHR.cov_ldsc_chr_wf as cov_ldsc_full{
            input:
                plink_format_prefix = plink_format_prefixes[index],
                cov_file = cov_file,
                out_prefix = output_file_prefixes[index],
                ldsc_cpu = ldsc_cpu,
                ldsc_mem_gb = ldsc_mem_gb,
        }
    }

}
