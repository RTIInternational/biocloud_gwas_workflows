import "biocloud_gwas_workflows/meta_analysis/metal/gwas/submodules/gwas_meta_utils.wdl" as UTILS

workflow preprocessing {
  Boolean comma_separated = false
  File gwas_results
  String study_basename
  Array[Int] chromosomes_to_keep

  Int variant_id_column
  Int chromosome_column
  Int pos_column
  Int noncoded_allele_column
  Int coded_allele_column
  Int effect_size_column
  Int standard_error_column
  Int pvalue_column
  String container_source

  String docker_ubuntu = if(container_source == "dockerhub") then "ubuntu:22.04-alpine" else "public.ecr.aws/ubuntu/ubuntu:22.04-alpine"
  String docker_python = if(container_source == "dockerhub") then "python:3.12-alpine" else "public.ecr.aws/docker/library/python:3.12-alpine"

  # Unzip file if it needs to be unzipped
  if(basename(gwas_results) != basename(gwas_results,".gz")){
    call UTILS.gunzip as gunzip {
      input:
        in_file = gwas_results,
        docker = docker_ubuntu
    }
    File unzipped_gwas_results = gunzip.output_file
  }

  # This small section is because we have to choose whether
  # to split the original input file or the unzipped input from previous task
  Array[File?] possible_files = [unzipped_gwas_results, gwas_results]
  File to_split = select_first(possible_files)

  # split the GWAS results up by chromosome
  call UTILS.split_by_chromosome as split_gwas {
    input:
      comma_separated = comma_separated,
      gwas_results = to_split,
      chromosome_column = chromosome_column,
      study_basename = study_basename,
      chromosomes_to_keep = chromosomes_to_keep,
      docker = docker_ubuntu
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
        pos_column = pos_column,
        noncoded_allele_column = noncoded_allele_column,
        coded_allele_column = coded_allele_column,
        effect_size_column = effect_size_column,
        standard_error_column = standard_error_column,
        pvalue_column = pvalue_column,
        docker = docker_ubuntu
    }
  }

  output {
    Array[File] proc_output = keep_columns.metal_input
    Array[Int] chromosome_order = split_gwas.chr_order
  }
}
