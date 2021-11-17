import "biocloud_gwas_workflows/meta_analysis/metal/ewas/submodules/ewas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/submodules/ewas_meta_preprocessing.wdl" as PREPROCESS 
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/submodules/ewas_meta_postprocessing.wdl" as POSTPROCESS 


workflow metal_ewas_meta_analysis_wf {
  Array[File] ewas_results
  Array[String] study_basename
  Array[Int] chromosomes_to_keep = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22] # default to autosomes
  String ancestry
  String plot_basename
  Float pvalue_threshold
  Boolean remove_singletons
  Float manhattan_significance_value = 6.3 # for ewas (log 1.22*10^{-7})

  Array[Int] probe_id_column
  Array[Int] position_column
  Array[Int] chromosome_column
  Array[Int] effect_size_column
  Array[Int] standard_error_column
  Array[Int] pvalue_column


  # Prepare input files for meta-analysis
  scatter (ewas_index in range(length(ewas_results))) {
    call PREPROCESS.preprocessing as preprocessing {
      input:
        ewas_results = ewas_results[ewas_index],
        study_basename = study_basename[ewas_index],
        chromosome_column = chromosome_column[ewas_index],
        chromosomes_to_keep = chromosomes_to_keep,

        probe_id_column = probe_id_column[ewas_index],
        pos_column = position_column[ewas_index],
        chromosome_column = chromosome_column[ewas_index],
        effect_size_column = effect_size_column[ewas_index],
        pvalue_column = pvalue_column[ewas_index],
        standard_error_column = standard_error_column[ewas_index],
    }
  }

  # transpose so we can group each set of chromosome files
  Array[Array[File]] metal_input = transpose(preprocessing.proc_output)

  # Perform meta-analysis using METAL on each chromosome
  scatter (chrom_order in range(length(metal_input))) {
    call UTILS.run_metal as metal {
      input:
        ewas_files = metal_input[chrom_order],
        chromosome = preprocessing.chromosome_order[0][chrom_order],
        ancestry = ancestry
    }
  }

  # Create figures and tables. 
  call POSTPROCESS.postprocessing as postprocessing {
    input:
      pvalue_threshold = pvalue_threshold,
      plot_basename = plot_basename,
      metal_results = metal.metal_results,
      remove_singletons = remove_singletons,
      manhattan_significance_value = manhattan_significance_value
  }

 output {
   File meta_analysis_full_results = postprocessing.finalTable
   File meta_analysis_top_hits = postprocessing.topHits
   Array[File] meta_analysis_plots = postprocessing.plots
  }
}

