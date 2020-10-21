import "biocloud_gwas_workflows/biocloud_wdl_tools/genesis/genesis.wdl" as GENESIS
# import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
# import "biocloud_gwas_workflows/biocloud_wdl_tools/make_gwas_summary_stats/make_gwas_summary_stats.wdl" as STAT
# import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
# import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
# import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

workflow genesis_gwas_chr_wf{
    File fileInGeno
    File fileInPheno
    File fileOut

    String genoFormat       # Options: gds
    String pheno            # Column name in phenotype file
    Array[String]? covars   # Array of column names of covars
    String family           # Options: gaussian
    String gxE              # Column name in phenotype file for GxE interaction
    Boolean? gzip = false

    # For annotating with population MAF info
    # File pop_maf_file
    # String maf_population

    Int genesis_cpu = 1
    Int genesis_mem_gb = 2

    call GENESIS.genesis{
        input:
            fileInGeno = fileInGeno,
            fileInPheno = fileInPheno,
            fileOut = fileOut,
            genoFormat = genoFormat,
            pheno = pheno,
            covars = covars,
            family = family,
            gxE = gxE,
            gzip = gzip,
            cpu = genesis_cpu,
            mem_gb = genesis_mem_gb
    }

    output {
        File summary_stats = fileOut
    }
}