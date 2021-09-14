import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS
import "/shared/jmarks/biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT


workflow postprocessing {
  String plot_basename
  Float pvalue_threshold
  Array[File] metal_results

  # exclude singletons
  scatter (chrom_results in range(length(metal_results))) {
    Int chr = sub(sub(metal_results[chrom_results], "^.+analysis_chr", ""), "_output1.metal", "")

    call UTILS.exclude_singletons as singletons {
      input:
        ewas_file = metal_results[chrom_results],
        chromosome = chr
    }
  }

  # merge the chromosome-specific results files into one chromosome-sorted file
  call UTILS.merge_results as merge {
    input:
      ewas_results = singletons.singletons_output
  }


  # create a P-value filtered table
  call UTILS.top_results as top_results {
    input:
      ewas_results = merge.merged_results,
      pvalue = pvalue_threshold
  }


  #Chromosome Position MarkerName Effect StdErr P-value Direction
  # plot resuls
  call PLOT.generate_gwas_plots as ewas_plots {
    input:
      summary_stats = merge.merged_results,
      col_id = "MarkerName",
      col_chromosome = "Chromosome",
      col_position = "Position",
      col_p = "P-value",
      output_basename = plot_basename
  }

  output {
   Array[File] finalTable = singletons.singletons_output
   File topHits = top_results.final_results_pvalue_filtered
   Array[File] plots = ewas_plots.plots
  }
}
