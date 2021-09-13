import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/ewas_meta_utils.wdl" as UTILS
import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/ewas_meta_preprocessing.wdl" as PREPROCESS 
import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/ewas_meta_postprocessing.wdl" as POSTPROCESS 


workflow metal_ewas_meta_analysis_wf {
  Array[File] ewas_results
  Array[String] study_basename
  Array[Int] chromosomes_to_keep
  String ancestry
  #String plot_basename

  Array[Int] probe_id_column
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
      #plot_basename = plot_basename,
      metal_results = metal.metal_results
  }
}
