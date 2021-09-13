import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/ewas_meta_utils.wdl" as UTILS

workflow preprocessing {
  File ewas_results
  String study_basename
  Array[Int] chromosomes_to_keep

  Int probe_id_column
  Int chromosome_column
  Int effect_size_column
  Int standard_error_column
  Int pvalue_column

  # split the EWAS results up by chromosome
  call UTILS.split_by_chromosome as split_ewas {
    input:
      ewas_results = ewas_results,
      study_basename = study_basename,
      chromosome_column = chromosome_column,
      chromosomes_to_keep = chromosomes_to_keep
  }
  
  Array[File] chr_split_results = split_ewas.chr_files
  Array[Int] kept_chroms = split_ewas.chr_order


  # keep only columns necessary for METAL 
  scatter (chrom_order in range(length(chr_split_results))) {
    call UTILS.keep_columns as keep_columns {
      input:
        infile =  chr_split_results[chrom_order],
        chromosome = kept_chroms[chrom_order],
        study_basename = split_ewas.ewas_name,

        probe_id_column = probe_id_column,
        chromosome_column = chromosome_column,
        effect_size_column = effect_size_column,
        pvalue_column = pvalue_column,
        standard_error_column = standard_error_column,
    }
  }

  output {
    Array[File] proc_output = keep_columns.metal_input
    Array[Int] chromosome_order = split_ewas.chr_order
  }
}
