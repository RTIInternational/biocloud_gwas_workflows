import "submodules/gwas_meta_utils.wdl" as UTILS
import "submodules/gwas_meta_preprocessing.wdl" as PREPROCESS 
import "submodules/gwas_meta_postprocessing.wdl" as POSTPROCESS 


workflow metal_gwas_meta_analysis_wf {

  Array[File] gwas_results = ["/home/ec2-user/rti-hiv/gwas_meta/hiv_acquisition/0034/data/mclaren_ea_hiv_acquisition_hg39_all_chr_rsq0.30_maf_0.01_rsid_only.tsv.gz", "/home/ec2-user/rti-hiv/gwas_meta/hiv_acquisition/0034/data/uhs1_ea_hiv_acquisition_all_chr_stats_rsq0.30_maf0.01.txt.gz"]
  Array[String] study_basename = ["mclaren", "uhs1"]
  Array[Int] chromosomes_to_keep  = [1, 2]
  String ancestry = "eur"
  String plot_basename = "jess-test"
  String remove_singletons = "true"

  Array[Int] variant_id_column = [1, 3]
  Array[Int] chromosome_column = [2, 1]
  Array[Int] pos_column = [3, 2]
  Array[Int] effect_size_column = [9, 12]
  Array[Int] coded_allele_column = [4, 4]
  Array[Int] noncoded_allele_column = [5, 5]
  Array[Int] pvalue_column = [11, 16]
  Array[Int] standard_error_column = [10, 13]
#  Array[Int] maf_column
#  Array[Int] rsquared_column


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
       standard_error_column = standard_error_column[gwas_index]
        #maf_column = maf_column[gwas_index],
        #rsquared_column = rsquared_column[gwas_index]  
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
        ancestry = ancestry

    }
  }

  # Create figures and tables. 
  call POSTPROCESS.postprocessing as postprocessing {
    input:
      metal_results = metal.metal_results,
      plot_basename = plot_basename,
      remove_singletons = remove_singletons

  }

  output {
   File meta_analysis_full_results = postprocessing.all_results
   File meta_analysis_top_hits = postprocessing.top_results
   Array[File] meta_analysis_plots = postprocessing.plots
  }
}
