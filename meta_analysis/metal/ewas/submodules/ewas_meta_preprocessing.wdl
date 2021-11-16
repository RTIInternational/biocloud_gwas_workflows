import "biocloud_gwas_workflows/meta_analysis/metal/ewas/submodules/ewas_meta_utils.wdl" as UTILS

workflow preprocessing {
  Boolean comma_separated = true
  File ewas_results
  String study_basename

  Int probe_id_column
  Int pos_column
  Int chromosome_column
  Int effect_size_column
  Int standard_error_column
  Int pvalue_column
  Array[Int] chromosomes_to_keep

  # Unzip file if it needs to be unzipped
  if(basename(ewas_results) != basename(ewas_results,".gz")){
    call UTILS.gunzip as gunzip {
      input:
        in_file = ewas_results
    }
    File unzipped_ewas_results = gunzip.output_file
  }

  # This small section is because we have to choose whether
  # to split the original input file or the unzipped input from previous task
  Array[File?] possible_files = [unzipped_ewas_results, ewas_results]
  File to_split = select_first(possible_files)

  # split the EWAS results up by chromosome
  call UTILS.split_by_chromosome as split_ewas {
    input:
      comma_separated = comma_separated,
      ewas_results = to_split,
      chromosome_column = chromosome_column,
      study_basename = study_basename,
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
        pos_column = pos_column,
        chromosome_column = chromosome_column,
        effect_size_column = effect_size_column,
        standard_error_column = standard_error_column,
        pvalue_column = pvalue_column,
    }
  }

  output {
    Array[File] proc_output = keep_columns.metal_input
    Array[Int] chromosome_order = split_ewas.chr_order
  }
}
