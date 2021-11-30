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
    Array[File] ewas_results = ["data/alspac_ea_ewas_results_chr21_2021-11-17.csv", "data/alspac_ea_ewas_results_chr22_2021-11-17.csv"]
    File prepare_table = "prepare_plot_table.R"
    Float fdr = 0.05

    String docker = "rtibiocloud/ewas:v0.0.1_fbfc0f1"

    command <<<

    Rscript ${prepare_table} \
        --input-files "${sep=" " ewas_results}" \
        --fdr ${fdr}

    >>>

    output {
        Float fdr_output = read_float("fdr_${fdr}_adjusted_gw_threshold.txt")
        Float bonferroni = read_float("bonferroni_adjusted_gw_threshold.txt")
        File plot_table = "plotting_table.csv"
    }

    runtime {
        docker: docker
        cpu: "1"
        memory: "4 GB"
    }
}

task plot_ewas {
    File plot_table = "/home/ec2-user/rti-cannabis/ewas/alspac/scripts/cromwell-executions/ewas_test/061f9466-b6c7-45aa-a761-b05f941b1ae7/call-plot_table/execution/plotting_table.csv"
    File plot_script = "make_plots.R"
    Array[String] colors = ["red", "blue"]
    String plot_basename = "jesseplot"
    Float bonferroni =   0.00000459
    Float fdr = 0.0000092
    String docker = "rtibiocloud/ewas:v0.0.1_fbfc0f1"
    

    command
    <<<

    Rscript ${plot_script} \
        --table "${plot_table}" \
        --fdr ${fdr} \
        --bonferroni ${bonferroni} \
        --colors "${sep=" " colors}" \
        --output "${plot_basename}"

    >>>

    output {
        Array[File] plots = glob("*png")
    }

    runtime {
        docker: docker
        cpu: "1"
        memory: "4 GB"
    }
}