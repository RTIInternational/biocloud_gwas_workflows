import "/home/ec2-user/wdl/biocloud_gwas_workflows/meta_analysis/metal/ewas/submodules/ewas_meta_utils.wdl" as UTILS
import "/home/ec2-user/wdl/biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT


workflow postprocessing {
  String plot_basename
  Float pvalue_threshold
  Boolean remove_singletons
  Float manhattan_significance_value # 7.3 for gwas and ~6.3 for ewas
  Array[File] metal_results
  #Array[File] metal_results = ["/home/ec2-user/wdl/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/cromwell-executions/metal_ewas_meta_analysis_wf/7e941f50-9c0b-4791-ac04-50ebaf6703dd/call-postprocessing/POSTPROCESS.postprocessing/76c481d3-69e3-4b3d-beb2-c6ef61bec428/call-prune/shard-0/inputs/-1239375637/afr_meta_analysis_chr1_output1.metal"]


  scatter (chrom_results in range(length(metal_results))) {
    Int chrom = sub(sub(metal_results[chrom_results], "^.+analysis_chr", ""), "_output1.metal", "")
    call UTILS.prune_results as prune {
      input:
        chromosome = chrom,
        ewas_file = metal_results[chrom_results]
    }
  }


  # exclude singletons
  if (remove_singletons) {
    scatter (chrom_results in range(length(metal_results))) {
      Int chr = sub(sub(metal_results[chrom_results], "^.+analysis_chr", ""), "_output1.metal", "")
      call UTILS.exclude_singletons as singletons {
        input:
          ewas_file = prune.pruned_columns[chrom_results],
          chromosome = chr
      }
    }
  }

  Array[File] ewas_results = select_first([singletons.singletons_output, prune.pruned_columns])

  # merge the chromosome-specific results files into one chromosome-sorted file
  call UTILS.merge_results as merge {
    input:
      ewas_results = ewas_results
  }


  # create a P-value filtered table
  call UTILS.top_results as top_results {
    input:
      ewas_results = merge.merged_results,
      pvalue = pvalue_threshold
  }


  # plot resuls
  call PLOT.generate_gwas_plots as ewas_plots {
    input:
      summary_stats = merge.merged_results,
      col_id = "MarkerName",
      col_chromosome = "Chromosome",
      col_position = "Position",
      col_p = "P-value",
      output_basename = plot_basename,
      manhattan_significance_value = manhattan_significance_value
  }

  output {
   #Array[File] finalTable = ewas_results
    File finalTable = merge.merged_results
    File topHits = top_results.final_results_pvalue_filtered
    Array[File] plots = ewas_plots.plots
  }
}
