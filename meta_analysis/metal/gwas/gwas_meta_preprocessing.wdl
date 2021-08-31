import "biocloud_gwas_workflows/meta_analysis/metal/gwas/gwas_meta_utils.wdl" as UTILS

workflow preprocessing {
  File gwas_results
  String study_basename
  Array[Int] chromosomes_to_keep

  Int variant_id_column
  Int chromosome_column
  Int effect_size_column
  Int coded_allele_column
  Int noncoded_allele_column
  Int pvalue_column
  Int standard_error_column
  Int maf_column
  Int rsquared_column


  # split the GWAS results up by chromosome
  call UTILS.split_by_chromosome as split_gwas {
    input:
      gwas_results = gwas_results,
      study_basename = study_basename,
      chromosome_column = chromosome_column,
      chromosomes_to_keep = chromosomes_to_keep
  }
  
  Array[File] chr_split_results = split_gwas.chr_files
  Array[Int] kept_chroms = split_gwas.chr_order


  # keep only columns necessary for METAL 
  scatter (chrom_order in range(length(chr_split_results))) {
    call UTILS.keep_columns as keep_columns {
      input:
        infile =  chr_split_results[chrom_order],
        chromosome = kept_chroms[chrom_order],
        study_basename = split_gwas.gwas_name,

        variant_id_column = variant_id_column,
        chromosome_column = chromosome_column,
        effect_size_column = effect_size_column,
        coded_allele_column = coded_allele_column,
        noncoded_allele_column = noncoded_allele_column,
        pvalue_column = pvalue_column,
        standard_error_column = standard_error_column,
        maf_column = maf_column,
        rsquared_column = rsquared_column,
    }
  }

  output {
    Array[File] proc_output = keep_columns.metal_input
    Array[Int] chromosome_order = split_gwas.chr_order
  }
}
