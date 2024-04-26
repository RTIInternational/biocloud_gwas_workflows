import "biocloud_gwas_workflows/meta_analysis/metal/gwas/submodules/gwas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/meta_analysis/metal/gwas/submodules/gwas_meta_preprocessing.wdl" as PREPROCESS 
import "biocloud_gwas_workflows/meta_analysis/metal/gwas/submodules/gwas_meta_postprocessing.wdl" as POSTPROCESS 


workflow metal_gwas_meta_analysis_wf {

  Array[File] gwas_results
  Array[String] study_basename
  Array[Int] chromosomes_to_keep
  String ancestry
  String full_results_name
  String remove_singletons
  Float pvalue_threshold = 0.001

  Array[Int] variant_id_column
  Array[Int] chromosome_column
  Array[Int] pos_column
  Array[Int] effect_size_column
  Array[Int] coded_allele_column
  Array[Int] noncoded_allele_column
  Array[Int] pvalue_column
  Array[Int] standard_error_column
#  Array[Int] maf_column
#  Array[Int] rsquared_column
  String container_source = "ecr" # or "dockerhub"

  # Prepare input files for meta-analysis
  scatter (gwas_index in range(length(gwas_results))) {
    call PREPROCESS.preprocessing as preprocessing {
      input:
        gwas_results = gwas_results[gwas_index],
        study_basename = study_basename[gwas_index],
        chromosome_column = chromosome_column[gwas_index],
        chromosomes_to_keep = chromosomes_to_keep,

        variant_id_column = variant_id_column[gwas_index],
        chromosome_column = chromosome_column[gwas_index],
        pos_column = pos_column[gwas_index],
        effect_size_column = effect_size_column[gwas_index],
        coded_allele_column = coded_allele_column[gwas_index],
        noncoded_allele_column = noncoded_allele_column[gwas_index],
        pvalue_column = pvalue_column[gwas_index],
        standard_error_column = standard_error_column[gwas_index],
        #maf_column = maf_column[gwas_index],
        #rsquared_column = rsquared_column[gwas_index]
        container_source = container_source
    }
  }

  # transpose so we can group each set of chromosome files
  Array[Array[File]] metal_input = transpose(preprocessing.proc_output)

  # Perform meta-analysis using METAL on each chromosome
  scatter (chrom_order in range(length(metal_input))){
    call UTILS.run_metal as metal {
      input:
        gwas_files = metal_input[chrom_order],
        chromosome = preprocessing.chromosome_order[0][chrom_order],
        ancestry = ancestry,
        container_source = container_source
    }
  }

  # Create figures and tables. 
  call POSTPROCESS.postprocessing as postprocessing {
    input:
      metal_results = metal.metal_results,
      full_results_name = full_results_name,
      remove_singletons = remove_singletons,
      pvalue_threshold = pvalue_threshold,
      container_source = container_source
  }

  output {
   File meta_analysis_full_results = postprocessing.all_results
   File meta_analysis_top_hits = postprocessing.top_results
   Array[File] meta_analysis_plots = postprocessing.plots
  }
}
