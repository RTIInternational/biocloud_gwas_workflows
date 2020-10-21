import "biocloud_gwas_workflows/association_testing/genesis/genesis_gwas_chr_wf.wdl" as GENESIS_CHR

workflow genesis_gwas_wf{
    Array[File] genotype_files
    Array[File] output_files
    File pheno_file
    String pheno_name
    # File? covar_file
    Array[String]? covars
    String geno_format
    String family
    String gxE
    Boolean? gzip = false

    Int genesis_cpu = 1
    Int genesis_mem_gb = 2

    # Do genesis chr workflow on each chromosome in parallel
    scatter(index in range(length(genotype_files))){

        call GENESIS_CHR.genesis_gwas_chr_wf as genesis{
            input:
                fileInGeno = genotype_files[index],
                fileOut = output_files[index],
                fileInPheno = pheno_file,
                genoFormat = geno_format,
                pheno = pheno_name,
                covars = covars,
                family = family,
                gxE = gxE,
                gzip = gzip,
                genesis_cpu = genesis_cpu,
                genesis_mem_gb = genesis_mem_gb
        }

    }
}