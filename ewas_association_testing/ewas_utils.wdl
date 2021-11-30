task run_ewas_rscript {

    File pheno_file
    File dnam_file
    File ewas_rscript
    String sample_name
    String output_basename

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

    parameter_meta  {
        pheno_file: "Phenotype file that contains the outcome as well as the covariates."
        dnam_file: "DNA methylation data post-QC. The sample names should match that given in the phenotype file."
        sample_name: "Sample name given in the DNAm data as well as the phenotype file."
        ewas_rscript: "Rscript that performs the EWAS."
        output_basename: "A descriptive basename for the output file--e.g. alspac_ea_ewas_model_1"
        docker: "Docker image to use for the EWAS. e.g. rtibiocloud/ewas:v0.0.1_fbfc0f1"
    }

    meta {
        description: "Performs an epigenome wide analysis EWAS using DNAm data and a phenotype file."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task plot_table {
    Array[File] ewas_results
    File prepare_table
    Float fdr

    String docker

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

    parameter_meta  {
        ewas_results: "An array of chromosome level EWAS results."
        prepare_table: "An Rscript that takes the chromosome level EWAS results and consolidates them into a single table for plotting."
        fdr: "False discovery rate."
        docker: "Docker image to use for the EWAS plotting table generation. e.g. rtibiocloud/ewas:v0.0.1_fbfc0f1"
    }

    meta {
        description: "Prepares the results from the EWAS into a table to create Manhattan and QQ plots."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task plot_ewas {
    File plot_table
    File plot_script
    Array[String] colors
    String plot_basename
    Float bonferroni
    Float fdr
    String docker

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

    parameter_meta  {

        plot_table: "Table for plotting. Contains p-values."
        plot_script: "Rscript that performs the Manhattan and QQ plotting."
        colors: "An Array of two colors for the Manhattan plot."
        plot_basename: "Basename for the Manhattan and QQ plots."
        bonferroni: "Bonferroni adjusted genome-wide significance threshold."
        fdr: "FDR adjusted genome-wide significance threshold (usually smaller than bonferroni). FDR of <some-specified threshold> is reach at this this value."
        docker: "Docker image to use for the EWAS plotting. e.g. rtibiocloud/ewas:v0.0.1_fbfc0f1"
    }

    meta {
        description: "Creates Manhattan and QQ plots."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}
