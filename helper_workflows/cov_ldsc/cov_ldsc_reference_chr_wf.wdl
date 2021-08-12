import "biocloud_gwas_workflows/biocloud_wdl_tools/cov_ldsc/cov_ldsc.wdl" as COVLDSC

workflow cov_ldsc_chr_wf{
    String plink_format_prefix
    File cov_file
    String out_prefix
    
    # Run cov-LDSC on the input file
    call COVLDSC.cov_ldsc as cov_ldsc{
        input:
            bfile = plink_format_prefix,
            cov_eigenvec = cov_file,
            out_prefix = out_prefix
            mem = ldsc_mem_gb
    }

    output {
        File mFile = cov_ldsc.m_File
        File m5_50_file = cov_ldsc.m_5_50_File
        File ldscore_out = cov_ldsc.ldscore_out
        File ldsc_log = cov_ldsc.logFile
    }
}

