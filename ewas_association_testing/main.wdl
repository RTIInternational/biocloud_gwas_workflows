import "biocloud_gwas_workflows/ewas_association_testing/single_ewas_wf.wdl" as EWAS
import "biocloud_gwas_workflows/ewas_association_testing/ewas_utils.wdl" as UTILS

workflow ewas_full {
    File pheno_file
    Array[File] dnam_files
    String test_var
    String sample_name
    Array[String] covariates
    String output_basename
    Array[String] plot_colors = ["red", "blue"] 
    Float fdr_value = 0.05

    String docker = "rtibiocloud/ewas:v0.0.2"

    scatter (dnam in dnam_files) {
        call EWAS.single_ewas as one_ewas {
            input:
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
            fdr = fdr_value,
            docker = docker
    }

    call UTILS.plot_ewas as plot {
        input:
            plot_table = table.plot_table,
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
