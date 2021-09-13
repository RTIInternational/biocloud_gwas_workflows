import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/ewas_meta_utils.wdl" as UTILS
#import "/shared/jmarks/biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT


workflow postprocessing {
#  String plot_basename = "plotty_plot"
#  Float pvalue_threshold = 0.001
#  Array[File] metal_results = [[""], [""]]
  Array[File] metal_results = ["/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/cromwell-executions/metal_ewas_meta_analysis_wf/224bd07e-f7be-4e50-8fe4-84b1e65bf465/call-metal/shard-0/execution/afr_meta_analysis_chr1_output1.metal"]

  # exclude singletons and add Chromosome & Position columns
  scatter (chrom_results in range(length(metal_results))) {
    Int chr = sub(sub(metal_results[chrom_results], "^.+analysis_chr", ""), "_output1.metal", "")

    call UTILS.exclude_singletons as singletons {
      input:
        ewas_file = metal_results[chrom_results],
        chromosome = chr
    }
  }


#  # merge the chromosome-specific results files into one chromosome-sorted file
#  call UTILS.merge_results as merge {
#    input:
#      ewas_results = singletons.singletons_output
#  }
#
#
#  # create a P-value filtered table
#  call UTILS.final_results as final_table {
#    input:
#      ewas_results = merge.merged_results,
#      pvalue = pvalue_threshold
#  }
#
#
#  # plot resuls
#  call PLOT.generate_gwas_plots as ewas_plots {
#    input:
#      summary_stats = merge.merged_results,
#      col_id = "MarkerName",
#      col_chromosome = "Chromosome",
#      col_position = "Position",
#      col_p = "P-value",
#      output_basename = plot_basename
#  }
#
  output {
   Array[File] remove_singletons_temp_return = singletons.singletons_output
   #File finalTable = final_table.final_table
   #Array[File] plots = ewas_plots.plots
  }
}
