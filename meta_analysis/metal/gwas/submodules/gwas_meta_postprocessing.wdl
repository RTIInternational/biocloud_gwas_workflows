import "biocloud_gwas_workflows/meta_analysis/metal/gwas/submodules/gwas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT


workflow postprocessing {
  String plot_basename 
  Float pvalue_threshold = 0.001
  Array[File] metal_results = ["/home/ec2-user/rti-hiv/gwas_meta/hiv_acquisition/0034/biocloud_gwas_workflows/meta_analysis/metal/gwas/cromwell-executions/metal_gwas_meta_analysis_wf/ab6f27c9-4d2a-45fb-82dc-dfcc3939750d/call-metal/shard-0/execution/eur_meta_analysis_chr1_output1.metal", "/home/ec2-user/rti-hiv/gwas_meta/hiv_acquisition/0034/biocloud_gwas_workflows/meta_analysis/metal/gwas/cromwell-executions/metal_gwas_meta_analysis_wf/ab6f27c9-4d2a-45fb-82dc-dfcc3939750d/call-metal/shard-1/execution/eur_meta_analysis_chr2_output1.metal"]
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
  call UTILS.merge_results as merge {
    input:
      gwas_results = singletons.singletons_output
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
      output_basename = plot_basename
  }

  output {
   File all_results = merge.merged_results
   File top_results = final_table.final_results_pvalue_filtered
   Array[File] plots = gwas_plots.plots
  }
}
