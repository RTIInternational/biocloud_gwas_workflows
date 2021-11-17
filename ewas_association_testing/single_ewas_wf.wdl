import "biocloud_gwas_workflows/ewas_association_testing/ewas_utils.wdl" as UTILS

workflow single_ewas {
    File pheno_file = "pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    File dnam_file = "alspac_dnameth_betas_chr21.rda"
    File ewas_rscript = "cannabis_alspac_ea_model1.R"
    String sample_name = "Sample_Name"
    String output_basename = "alspac_ewas_results" # e.g., alspac_ewas_results ==> alspac_ewas_results_chr21.csv

    String docker = "ffang8/ewas:v041221"

    call UTILS.run_ewas_rscript as rscript {
        input:
            pheno_file = pheno_file,
            dnam_file = dnam_file,
            ewas_rscript = ewas_rscript,
            sample_name = sample_name,
            output_basename = output_basename,
            
            docker = docker
    }

    output {
        Array[File] ewas_output = rscript.ewas_results # this is actually just a single file but because of glob...
    }
}