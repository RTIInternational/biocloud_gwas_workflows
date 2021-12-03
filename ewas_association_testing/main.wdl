import "biocloud_gwas_workflows/ewas_association_testing/single_ewas_wf.wdl" as EWAS
import "biocloud_gwas_workflows/ewas_association_testing/ewas_utils.wdl" as UTILS

workflow ewas_full {
    File ewas_rscript #= "cannabis_alspac_ea_model1.R"
    File pheno_file #= "data/pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    Array[File] dnam_files #= ["data/alspac_dnameth_betas_chr21.rda"]
    String test_var #= "cannabisUse"
    String sample_name #= "Sample_Name"
    Array[String] covariates #= ["age_at_DNAm", "sv1", "sv2", "sv3", "sv4", "sv5", "sv6", "Bcell", "CD4T", "CD8T", "Mono", "Gran", "NK"]
    String output_basename #= "alspac_ea_ewas_results" ==> alspac_ewas_results_chr21.csv
    Array[String] plot_colors #= ["red", "blue"] 
    Float fdr_value #= 0.05
    File plot_script #= "make_plots.R"
    File prepare_table #= "prepare_plot_table.R"

    String docker #= "rtibiocloud/ewas:v0.0.1_fbfc0f1"

    scatter (dnam in dnam_files) {
        call EWAS.single_ewas as one_ewas {
            input:
                ewas_rscript = ewas_rscript,
                pheno_file = pheno_file,
                dnam_file = dnam,
                test_var = test_var,
                sample_name = sample_name,
                covariates = covariates,
                output_basename = output_basename,
                docker = docker
        }
    }

    call UTILS.plot_table as table {
        input:
            ewas_results = one_ewas.ewas_output,
            prepare_table = prepare_table,
            fdr = fdr_value,
            docker = docker
    }

    call UTILS.plot_ewas as plot {
        input:
            plot_table = table.plot_table,
            plot_script = plot_script,
            colors = plot_colors,
            plot_basename = output_basename,
            fdr = table.fdr_output,
            bonferroni = table.bonferroni,
            docker = docker
    }

    output {
        Array[File] final_output = one_ewas.ewas_output
        Array[File] final_plots = plot.plots
    }
}
