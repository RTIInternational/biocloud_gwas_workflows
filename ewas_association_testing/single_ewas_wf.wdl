import "biocloud_gwas_workflows/ewas_association_testing/ewas_utils.wdl" as UTILS

workflow single_ewas {
    File pheno_file
    File dnam_file
    String test_var
    String sample_name
    Array[String] covariates
    String output_basename
    String docker = "rtibiocloud/ewas:v0.0.2"

    call UTILS.run_ewas_rscript as rscript {
        input:
            pheno_file = pheno_file,
            dnam_file = dnam_file,
            test_var = test_var,
            sample_name = sample_name,
            covariates = covariates,
            output_basename = output_basename,
            docker = docker
    }

    output {
        File ewas_output = rscript.ewas_results[0] # this is actually just a single file but because of glob it's an array. Ergo return 1st element of array.
    }
}
