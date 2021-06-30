import "biocloud_gwas_workflows/meta_analysis/metal/gwas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/meta_analysis/metal/gwas_meta_preprocessing.wdl" as PREPROCESS 
import "biocloud_gwas_workflows/meta_analysis/metal/gwas_meta_postprocessing.wdl" as POSTPROCESS 


workflow metal_gwas_meta_analysis_wf {
  Array[File] gwas_results
  Array[String] study_basename
  Array[Int] chromosomes_to_keep
  String ancestry
  String plot_basename

  Array[Int] variant_id_column
  Array[Int] chromosome_column
  Array[Int] effect_size_column
  Array[Int] coded_allele_column
  Array[Int] noncoded_allele_column
  Array[Int] pvalue_column
  Array[Int] standard_error_column
  Array[Int] maf_column
  Array[Int] rsquared_column



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
        effect_size_column = effect_size_column[gwas_index],
        coded_allele_column = coded_allele_column[gwas_index],
        noncoded_allele_column = noncoded_allele_column[gwas_index],
        pvalue_column = pvalue_column[gwas_index],
        standard_error_column = standard_error_column[gwas_index],
        maf_column = maf_column[gwas_index],
        rsquared_column = rsquared_column[gwas_index]  
    }
  }

  # transpose so we can group each set of chromosome files
  Array[Array[File]] metal_input = transpose(preprocessing.proc_output)

  # Perform meta-analysis using METAL on each chromosome
  scatter (chrom_order in range(length(metal_input))){

    #Int chr = sub(sub(metal_input[chrom_order][0], "^.+_chr", ""), "_specific_columns.tsv", "")

    call UTILS.run_metal as metal {
      input:
        gwas_files = metal_input[chrom_order],
        #chromosome = chr,
        chromosome = preprocessing.chromosome_order[0][chrom_order],
        ancestry = ancestry

    }
  }

  # Create figures and tables. 
  call POSTPROCESS.postprocessing as postprocessing {
    input:
      metal_results = metal.metal_results,
      plot_basename = plot_basename
  }
}
