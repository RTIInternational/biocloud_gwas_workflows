import "biocloud_gwas_workflows/biocloud_wdl_tools/genesis/genesis.wdl" as GENESIS
# import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
# import "biocloud_gwas_workflows/biocloud_wdl_tools/make_gwas_summary_stats/make_gwas_summary_stats.wdl" as STAT
# import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
# import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
# import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

workflow genesis_gwas_chr_wf{
    File fileInGeno
    File fileInPheno
    String outBaseName

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

    Boolean split_genotypes = true

    # Optionally chunk genotype file for parallel processing
    if(split_genotypes){
        # call SPLIT.split_vcf as split_vcf{
        #     input:
        #         input_vcf = vcf_in,
        #         records_per_split = records_per_split,
        #         output_basename = output_basename,
        #         cpu = split_vcf_cpu
        # }
    }

    # Loop through splits and do association testing on each
    # Array[File] split_genos = select_first([split_vcf.split_genos, [vcf_in]])
    Array[File] split_genos = [fileInGeno]

    scatter(split_index in range(length(split_genos))){

        # String split_output_basename = basename(sub(split_genos[split_index], "\\.gz$",""), ".vcf")
        String fileOut = outBaseName + ".assoc.tsv"

        # Run genesis for association
        call GENESIS.genesis{
            input:
                fileInGeno = split_genos[split_index],
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

        # # Remove header from association output file
        # File assoc_file = genesis.assoc_outputs[0]
        # call strip_genesis_headers{
        #     input:
        #         sumstats_in = assoc_file,
        #         output_basename = basename(assoc_file, ".gz")
        # }
    }

    output{
        File summary_stats = fileOut
    }

    # # Gather chunked genesis output if > 1 split
    # if(length(split_genos) > 1){
    #     # Collect chunked sumstats files into single zip folder
    #     call COLLECT.collect_large_file_list_wf as collect_genesis_sumstats{
    #         input:
    #             input_files = strip_genesis_headers.sumstats_out,
    #             output_dir_name = output_basename + ".genesis_output"
    #     }

    #     # Concat all sumstats files into single sumstat file
    #     call TSV.tsv_append as cat_genesis_sumstats{
    #         input:
    #             tsv_inputs_tarball = collect_genesis_sumstats.output_dir,
    #             output_filename = output_basename + ".genesis.MetaAssoc.tsv"
    #     }
    # }

    # # Use the combined sumstats file if >1 splits was merged, otherwise just use output from strip_headers call
    # File full_sumstats = select_first([cat_genesis_sumstats.tsv_output, strip_genesis_headers.sumstats_out[0]])

    # # Annotate sumstats with features from info file (R2, MAF) and pop MAF file (MAF from pop of interest)
    # call STAT.make_gwas_summary_stats as annotate_sumstats{
    #     input:
    #         file_in_summary_stats = full_sumstats,
    #         file_in_info = info_in,
    #         file_in_pop_mafs = pop_maf_file,
    #         file_in_summary_stats_format = "genesis",
    #         population = maf_population,
    #         file_out_prefix = output_basename + ".formatted"
    # }

    # output{
    #     File summary_stats = annotate_sumstats.output_file
    # }
}