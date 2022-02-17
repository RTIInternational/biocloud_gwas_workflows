import "biocloud_gwas_workflows/meta_analysis/metal/gwas/submodules/gwas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT


workflow postprocessing {
  String full_results_name
  Float pvalue_threshold = 0.001
  Array[File] metal_results
  String remove_singletons

  # exclude singletons and capitalize alleles
  scatter (chrom_results in range(length(metal_results))) {
    Int chr = sub(sub(metal_results[chrom_results], "^.+sis_chr", ""), "_output1.metal", "")
    call UTILS.exclude_singletons as singletons {
      input:
        gwas_file = metal_results[chrom_results],
        chromosome = chr,
        remove_singletons = remove_singletons
    }
  }


  # merge the chromosome-specific results files into one chromosome-sorted file
  call UTILS.merge_final_results as merge {
    input:
      gwas_results = singletons.singletons_output,
      full_results_name = full_results_name
  }


#  # create a P-value filtered table
  call UTILS.top_results as final_table {
    input:
      gwas_results = merge.merged_results,
      pvalue = pvalue_threshold
  }


  # plot resuls
  call PLOT.generate_gwas_plots as gwas_plots {
    input:
      summary_stats = merge.merged_results,
      col_id = "MarkerName",
      col_chromosome = "Chromosome",
      col_position = "Position",
      col_p = "P-value",
      output_basename = full_results_name
  }

  output {
   File all_results = merge.merged_results
   File top_results = final_table.final_results_pvalue_filtered
   Array[File] plots = gwas_plots.plots
  }
}
