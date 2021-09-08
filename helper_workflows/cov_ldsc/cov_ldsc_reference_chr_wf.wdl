import "biocloud_gwas_workflows/biocloud_wdl_tools/cov_ldsc/cov_ldsc.wdl" as COVLDSC

workflow cov_ldsc_chr_wf{
    File bed_in
    File bim_in
    File fam_in
    File cov_file
    String out_prefix
    Int covldsc_cpu = 16
    Int covldsc_mem_gb = 128

    # Run cov-LDSC on the input file
    call COVLDSC.cov_ldsc as cov_ldsc{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            cov_eigenvec = cov_file,
            out_prefix = out_prefix,
            cpu = covldsc_cpu,
            mem_gb = covldsc_mem_gb
    }

    output {
        File mFile = cov_ldsc.m_File
        File m5_50_file = cov_ldsc.m_5_50_File
        File ldscore_out = cov_ldsc.ldscore_out
        File ldsc_log = cov_ldsc.logFile
    }
}

