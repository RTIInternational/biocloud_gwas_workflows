import "biocloud_gwas_workflows/ewas_association_testing/single_ewas_wf.wdl" as EWAS
import "biocloud_gwas_workflows/ewas_association_testing/ewas_utils.wdl" as UTILS

workflow ewas_full {
    File pheno_file #= "data/pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt"
    Array[File] dnam_files #= ["data/alspac_dnameth_betas_chr21.rda"]
    String sample_name #= "Sample_Name"
    String output_basename #= "alspac_ea_ewas_results" ==> alspac_ewas_results_chr21.csv
    Array[String] plot_colors #= ["red", "blue"] 
    Float fdr_value #= 0.05
    File plot_script
    File ewas_rscript #= "cannabis_alspac_ea_model1.R"

    String docker = "rtibiocloud/ewas:v0.0.1_fbfc0f1"

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

    call UTILS.plot_table as plot_table {
        input:
            ewas_results = one_ewas.ewas_output,
            plot_script = plot_script,
            fdr_value = fdr_value,
            colors = plot_colors,
            plot_basename = output_basename,

            docker = "rtibiocloud/ewas:v0.0.1_fbfc0f1"
            #plot_script = "prepare_plots.R",
            #ewas_results = ["/home/ec2-user/rti-cannabis/ewas/alspac/scripts/cromwell-executions/ewas_full/9083f770-339f-42ce-971a-213f1f1fcb30/call-one_ewas/shard-0/EWAS.single_ewas/b28c14f3-55ee-4513-9ee4-f5feef1df838/call-rscript/execution/glob-55f68f231ff6123133f73be7360d3d4c/alspac_ea_ewas_results_chr21_2021-11-30.csv", "/home/ec2-user/rti-cannabis/ewas/alspac/scripts/cromwell-executions/ewas_full/9083f770-339f-42ce-971a-213f1f1fcb30/call-one_ewas/shard-1/EWAS.single_ewas/40edd8dc-487a-44aa-af5c-6157f0ee828c/call-rscript/execution/glob-55f68f231ff6123133f73be7360d3d4c/alspac_ea_ewas_results_chr22_2021-11-30.csv"],
    }

    output {
        Array[File] final_output = one_ewas.ewas_output
        Array[File] final_plots = plot_table.output_plots
    }
}
