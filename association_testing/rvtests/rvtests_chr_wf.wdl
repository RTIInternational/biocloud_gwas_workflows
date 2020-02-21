import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/rvtests/rvtests.wdl" as RV
import "biocloud_gwas_workflows/biocloud_wdl_tools/make_gwas_summary_stats/make_gwas_summary_stats.wdl" as STAT
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

workflow rvtests_chr_wf{
    File vcf_in
    File info_in
    File pheno_file
    String pheno_name
    File covar_file
    Array[String] covars
    File? kinship_matrix

    # For annotating with population MAF info
    File pop_maf_file
    String maf_population

    # Optional for X chr
    File? xHemiKinship_matrix
    Boolean? xHemi

    # Number of records per VCF split
    Int records_per_split = 50000

    # Set output basename
    String? user_output_basename
    String default_output_basename = basename(sub(vcf_in, "\\.gz$",""), ".vcf")
    String output_basename = select_first([user_output_basename, default_output_basename])

    # Chunk VCF and INFO files for parallel processing
    call SPLIT.split_vcf as split_vcf{
        input:
            input_vcf = vcf_in,
            records_per_split = records_per_split,
            output_basename = output_basename
    }

    # Loop through splits and do association testing on each
    scatter(split_index in range(length(split_vcf))){

        # Run rvtests for association
        call RV.rvtests{
            input:
                inVCF = split_vcf[split_index],
                phenoFile = pheno_file,
                out = output_basename,
                covarFile = covar_file,
                phenoName = pheno_name,
                covarsMaybe = covars,
                kinship = kinship_matrix,
                xHemiKinship = xHemiKinship_matrix,
                xHemi = xHemi
        }
    }

    # Collapse association outputs into single array (each RVTests produces an array of association files)
    # Pretty much just turns a 2-D file array into a 1-D file array
    call UTILS.flatten_string_array as flatten_sumstats{
        input:
            array = rvtests.assoc_outputs
    }

    # Collect chunked sumstats files into single zip folder
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
            input_files = rvtests.assoc_outputs,
            output_dir_name = output_basename + "rvtests_output"
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = output_basename + ".rvtests.merged_assoc.tsv"
    }

    call STAT.make_gwas_summary_stats as annotate_sumstats{
        input:
            file_in_summary_stats = cat_sumstats.tsv_output,
            file_in_info = info_in,
            file_in_pop_mafs = pop_maf_file,
            file_in_summary_stats_format = "rvtests",
            population = maf_population,
            file_out_prefix = output_basename + ".formatted"
    }

    output{
        File summary_stats = annotate_sumstats.output_file
    }
}