import "ewas_utils.wdl" as UTILS

workflow single_ewas {
    File pheno_file = "data/pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    File dnam_file = "data/alspac_dnameth_betas_chr21.rda"
    File ewas_rscript = "cannabis_alspac_ea_model1.R"
    String sample_name = "Sample_Name"
    String output_basename = "alspac_ewas_results" # e.g., alspac_ewas_results ==> alspac_ewas_results_chr21.csv

    String docker = "rtibiocloud/ewas:v0.0.1_fbfc0f1"

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
        File ewas_output = rscript.ewas_results[0] # this is actually just a single file but because of glob it's an array. Ergo return 1st element of array.
    }
}
