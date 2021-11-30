task run_ewas_rscript {

    File pheno_file #= "pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    File dnam_file #= "alspac_dnameth_betas_chr21.rda"
    File ewas_rscript #= "cannabis_alspac_ea_model1.R"
    String sample_name #= "Sample_Name"
    String output_basename #= "alspac_ewas_results" # e.g., alspac_ewas_results ==> alspac_ewas_results_chr21.csv

    String docker

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

task plot_table {
    Array[File] ewas_results #= ["data/alspac_ea_ewas_results_chr21_2021-11-17.csv", "data/alspac_ea_ewas_results_chr22_2021-11-17.csv"]

    File plot_script = "prepare_plots.R"
    Float fdr_value = 0.05
    Array[String] colors = ["red", "blue"]
    String plot_basename = "jesseplot"

    String docker = "rtibiocloud/ewas:v0.0.1_fbfc0f1"

    command <<<

    # eventually split this up so that only the plot table is created
    Rscript ${plot_script} \
        --input-files "${sep=" " ewas_results}" \
        --fdr ${fdr_value} \
        --colors "${sep=" " colors}" \
        --output "${plot_basename}"

    >>>

    output {
        #File out_table = "plotting_table.txt"
        Array[File] output_plots = glob("*png")
    }

    runtime {
        docker: docker
        cpu: "1"
        memory: "4 GB"
    }
}

#task plot_ewas {
#    
#}
