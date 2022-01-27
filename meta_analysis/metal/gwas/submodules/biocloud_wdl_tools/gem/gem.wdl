task gem {

    # Input/Output File Options:
    File pheno_file
    String out
    File? bgen
    File? sample
    File? pgen
    File? pvar
    File? psam
    File? bed
    File? bim
    File? fam
    String? output_style
    String out_log = out + ".log"

    # Phenotype File Options:
    String sampleid_name
    String pheno_name
    String? exposure_names
    String? int_covar_names
    String? covar_names
    Int? robust
    Float? tol
    String? delim
    String? missing_value
    Int? center
    Int? scale
    String? categorical_names
    Int? cat_threshold

    # Filtering Options:
    Float? maf
    Float? miss_geno_cutoff
    String? include_snp_file

    # Performance Options:
    Int? threads
    Int? stream_snps

    # Runtime attributes
    String docker = "rtibiocloud/gem:v1.4.2_05922f6"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3


    command {
        GEM --pheno-file ${pheno_file} \
            --out ${out} \
            --sampleid-name ${sampleid_name} \
            --pheno-name ${pheno_name} \
            ${"--bgen " + bgen} \
            ${"--sample " + sample} \
            ${"--pgen " + pgen} \
            ${"--pvar " + pvar} \
            ${"--psam " + psam} \
            ${"--bed " + bed} \
            ${"--bim " + bim} \
            ${"--fam " + fam} \
            ${"--output-style " + output_style} \
            ${"--exposure-names " + exposure_names} \
            ${"--int-covar-names " + int_covar_names} \
            ${"--covar-names " + covar_names} \
            ${"--robust " + robust} \
            ${"--tol " + tol} \
            ${"--delim " + delim} \
            ${"--missing-value " + missing_value} \
            ${"--center " + center} \
            ${"--scale " + scale} \
            ${"--categorical-names " + categorical_names} \
            ${"--cat-threshold " + cat_threshold} \
            ${"--maf " + maf} \
            ${"--miss-geno-cutoff " + miss_geno_cutoff} \
            ${"--include-snp-file " + include_snp_file} \
            ${"--threads " + threads} \
            ${"--stream-snps " + stream_snps} > \
            ${out_log}
    }

    output {
        File summary_stats = "${out}"
        File log = "${out_log}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

}
