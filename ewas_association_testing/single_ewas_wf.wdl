workflow single_ewas {
    File pheno_file = "pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    File dnam_file = "alspac_dnameth_betas_chr21.rda"
    File ewas_rscript = "cannabis_alspac_ea_model1.R"
    String sample_name = "Sample_Name"
    String output_basename = "alspac_ewas_results" # e.g., alspac_ewas_results ==> alspac_ewas_results_chr21.csv

    String docker = "ffang8/ewas:v041221"

    call run_ewas_rscript {
        input:
            pheno_file = pheno_file,
            dnam_file = dnam_file,
            ewas_rscript = ewas_rscript,
            sample_name = sample_name,
            output_basename = output_basename,
            
            docker = docker
    }

}

task run_ewas_rscript {

    File pheno_file #= "pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    File dnam_file #= "alspac_dnameth_betas_chr21.rda"
    File ewas_rscript #= "cannabis_alspac_ea_model1.R"
    String sample_name #= "Sample_Name"
    String output_basename #= "alspac_ewas_results" # e.g., alspac_ewas_results ==> alspac_ewas_results_chr21.csv

    String docker #= "ffang8/ewas:v041221"

    command 
    <<<
        Rscript ${ewas_rscript} \
        --pheno ${pheno_file} \
        --sample_name ${sample_name} \
        --methRdata ${dnam_file} \
        --output ${output_basename} 
    >>>

    output {
       Array[File] ewas_results = glob("*csv")
    }

    runtime {
        docker: docker
        cpu: "1"
        memory: "1 GB"
    }
}
