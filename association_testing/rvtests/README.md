# Association Testing with Rvtests among Unrelated Individuals

## Introduction

This document details the standard analysis workflow for performing association testing for unrelated individuals in Rvtests for binary or linear outcomes. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes genotypes in vcf format and phenotypes and covariates in Rvtests format as inputs and produces the following outputs:

1. Summary statistics from the association testing.
2. Summary statistics from the association testing for all variants with p < N, where N is a p-value threshold set by the user.
3. Q-Q plot of the association testing p-values.
4. Manhattan plot  of the association testing p-values.

The input and output formats are fully described in the appendix of this document.

Optional minor allele frequency (MAF) and imputation quality (RSQ) filters can be specified to limit the reported variants.

The steps in this workflow are as follows:

1. Perform association testing in Rvtests.
2. Add chr, position, and standard error to Rvtests results and convert to standard format for GWAS results (see appendix).
3. (Optional) Convert variant IDs to a different format.
4. (Optional) Filter by MAF in reference population.
5. (Optional) Filter by MAF in study population.
6. (Optional) Filter by imputation RSQ.
7. Filter summary statistics by p-value.
8. Generate Q-Q and Manhattan plots.

Each of these steps in described in detail below. Steps 1-6 are run separately for each chromosome. Steps 7 and 8 merge all chromosomes together.


### 1. Perform association testing in Rvtests.

Sample command for association testing for binary trait:
``` shell
rvtest \
    --inVcf [INPUT_VCF_FILE] \
    --pheno [INPUT_PHENO_FILE] \
    --pheno-name [PHENO_COL] \
    --covar [INPUT_COVAR_FILE] \
    --covar-name [COVAR_COL_1],...,[COVAR_COL_X] \
    --meta score \
    --dosage DS \
    --out [OUTPUT_PREFIX]
```

Sample command for association testing for continuous trait:
``` shell
rvtest \
    --noweb \
    --inVcf [INPUT_VCF_FILE] \
    --pheno [INPUT_PHENO_FILE] \
    --pheno-name [PHENO_COL] \
    --covar [INPUT_COVAR_FILE] \
    --covar-name [COVAR_COL_1],...,[COVAR_COL_X] \
    --meta score \
    --dosage DS \
    --useResidualAsPhenotype \
    --inverseNormal \
    --qtl \
    --out [OUTPUT_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `INPUT_VCF_FILE` | Genotype file in standard VCF format (see appendix) |
| `INPUT_PHENO_FILE` | Phenotype file in standard Rvtest format (see appendix) |
| `INPUT_COVAR_FILE` | Covariate file in standard Rvtest format (see appendix); may be the same as PHENO file |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUTPUT_PREFIX].MetaScore.assoc.gz` | Association testing results |
| `[OUTPUT_PREFIX].MetaScore.assoc.gz.tbi` | tabix index file for association test results |
| `[OUTPUT_PREFIX].log` | Log file for rvtest run |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--inVcf [INPUT_VCF_FILE]` | Genotypes in VCF format (see appendix) |
| `--pheno [INPUT_PHENO_FILE]` | Phenotypes in rvtest format (see appendix) |
| `--covar [INPUT_COVAR_FILE]` | Covariates in rvtest format (see appendix) |
| `--out [OUTPUT_PREFIX]` | Prefix for output files |
| `--pheno-name [PHENO_COL]` | Name of the column in INPUT_PHENO_FILE to use as outcome for association testing |
| `--covar-name [COVAR_COL_1],...,[COVAR_COL_X]` | Name of the column(s) in INPUT_COVAR_FILE to use as covariates for association testing |
| `--meta [MODEL]` | Model for association testing (score, dominant, or recessive) |
| `--dosage [DOSAGE_TAG]` | Tag for dosage in VCF file |
| `--useResidualAsPhenotype` | Fit regression model of phenotype on covariates and use residuals for association testing; typically used for continuous traits; usually used with `--inverseNormal` parameter |
| `--inverseNormal` | Inverse normalize residuals generated with `--useResidualAsPhenotype`; typically used for continuous traits |
| `--qtl` | Parameter indicating that the phenotype is a continuous trait |s
| `--noweb` | Skip check for remote version |


Full documentation for Rvtests can be found [here](http://zhanxw.github.io/rvtests/).


### 2. Add ID, MAF and standard error to Rvtests results and convert to standard format for GWAS results (see appendix).

Need to make formal script for this step.

Sample command:
``` shell
```


### 3. (Optional) Convert variant IDs to a different format.

Sample command for converting IDs to IMPUTE2 format:
``` shell
convert_to_1000g_p3_ids.pl \
    --file_in [INPUT_FILE] \
    --legend [LEGEND_FILE] \
    --file_out [OUTPUT_FILE] \
    --file_in_header 1 \
    --file_in_id_col 0 \
    --file_in_chr_col 1 \
    --file_in_pos_col 2 \
    --file_in_a1_col 3 \
    --file_in_a2_col 4 \
    --chr 22
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| INPUT_FILE | File in which to convert IDs; Can be any file that has columns for ID, chr, position, allele 1 and allele 2 |
| LEGEND_FILE | Reference panel legend file with desired IDs |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| OUTPUT_FILE | File identical to the INPUT_FILE, except that the contents of the ID column have been updated |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--file_in_header [N]` | Number of header rows in INPUT_FILE |
| `--file_in_id_col [N]` | ID column in INPUT_FILE (numbering starts with 0) |
| `--file_in_chr_col [N]` | Chromosome column in INPUT_FILE (numbering starts with 0) |
| `--file_in_pos_col [N]` | Position column in INPUT_FILE (numbering starts with 0) |
| `--file_in_a1_col [N]` | Allele 1 column in INPUT_FILE (numbering starts with 0) |
| `--file_in_a2_col [N]` | Allele 2 column in INPUT_FILE (numbering starts with 0) |
| `--chr [N]` | Chromosome represented in INPUT_FILE |


### 4. (Optional) Filter by MAF in reference population.

Sample command:
``` shell
/shared/bioinformatics/software/perl/utilities/extract_rows.pl \
    --source [INPUT_FILE] \
    --id_list [REMOVE_LIST] \
    --out [OUTPUT_FILE] \
    --header 1 \
    --id_column 0 \
    --remove
```


### 5. (Optional) Filter by MAF in study population.

Need to make formal script for this step.

Sample command:
``` shell
```


### 6. (Optional) Filter by imputation RSQ.

Need to make formal script for this step.

Sample command:
``` shell
```


### 7. Filter summary statistics by p-value.

Sample command:
``` shell
```


### 8. Generate Q-Q and Manhattan plots.

Sample command:
``` shell
/shared/bioinformatics/software/R/generate_gwas_plots.R \
    --in [INPUT_FILE] \
    --out [OUTPUT_FILE] \
    --in_chromosomes autosomal \
    --in_header \
    --col_id VARIANT_ID \
    --col_chromosome CHR \
    --col_position POSITION \
    --col_p P \
    --col_variant_type TYPE \
    --generate_snp_indel_manhattan_plot \
    --manhattan_odd_chr_color red \
    --manhattan_even_chr_color blue \
    --manhattan_points_cex 1.5 \
    --generate_snp_indel_qq_plot \
    --qq_lines \
    --qq_points_bg black \
    --qq_lambda
```
