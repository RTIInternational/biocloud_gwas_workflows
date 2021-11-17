import "biocloud_gwas_workflows/ewas_association_testing/single_ewas_wf.wdl" as EWAS

workflow ewas_all {
    File pheno_file #= "pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    Array[File] dnam_files #= ["alspac_dnameth_betas_chr21.rda"]
    File ewas_rscript #= "cannabis_alspac_ea_model1.R"
    String sample_name #= "Sample_Name"
    String output_basename #= "alspac_ea_ewas_results" # e.g., alspac_ewas_results ==> alspac_ewas_results_chr21.csv

    String docker #= "ffang8/ewas:v041221"

    scatter (dnam in dnam_files) {
        call EWAS.single_ewas as one_ewas {
            input:
                pheno_file = pheno_file,
                dnam_file = dnam,
                ewas_rscript = ewas_rscript,
                sample_name = sample_name,
                output_basename = output_basename,
                
                docker = docker
        }
    }

    output {
        Array[Array[File]] final_output = one_ewas.ewas_output
    }
}