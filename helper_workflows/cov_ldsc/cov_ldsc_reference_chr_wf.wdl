import "biocloud_gwas_workflows/biocloud_wdl_tools/cov_ldsc/cov_ldsc.wdl" as COVLDSC

workflow cov_ldsc_chr_wf{
    String plink_format_prefix
    File cov_file
    String out_prefix
    Int ldsc_cpu
    Int ldsc_mem_gb
    
    # Do cov-LDSC on the input file
    call COVLDSC.cov_ldsc as cov_ldsc{
        input:
            plink_format_prefix = plink_format_prefix,
            cov_file = cov_file,
            out_prefix = out_prefix
            cpu = ldsc_cpu,
            mem_gb = ldsc_mem_gb
    }

    output {
        File ldsc_out1 = cov_ldsc.out1
        File ldsc_out2 = cov_ldsc.out2
        File ldsc_out3 = cov_ldsc.out3
        File ldsc_log = cov_ldsc.log
    }
}

