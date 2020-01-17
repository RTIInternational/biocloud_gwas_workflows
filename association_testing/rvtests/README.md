# Association Testing with Rvtests among Unrelated Individuals

## Introduction

This living document details the Center for Omics Discovery and Epidemiology (CODE) analysis workflow for performing association testing for unrelated individuals in Rvtests for binary or linear outcomes. It is intended to serve as the standard data processing workflow and reference guide across CODE projects. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes genotypes in vcf format and phenotypes and covariates in Rvtests format as inputs and produces the following outputs:

1. Summary statistics from the association testing.
2. Summary statistics from the association testing for all variants with p < N, where N is a p-value threshold set by the user.
3. Q-Q plot of the association testing p-values.
4. Manhattan plot  of the association testing p-values.

The input and output formats are fully described in the appendix of this document.

Optional minor allele frequency (MAF) and imputation quality (RSQ) filters can be specified to limit the reported variants.

The steps in this workflow are as follows:

1. Perform association testing in Rvtests.
2. Add chr, position, and standard error to Rvtests results.
3. (Optional) Convert variant IDs to a different format.
4. (Optional) Filter by MAF in reference population.
5. (Optional) Filter by MAF in study population.
6. (Optional) Filter by imputation RSQ.
7. Generate Q-Q and Manhattan plots.
8. Filter summary statistics by p-value.

Each of these steps in described in detail below.


### 1. Perform association testing in Rvtests.

Sample command for association testing for binary trait:
``` shell
rvtest \
    --inVcf [INPUT_VCF_FILE] \
    --pheno [INPUT_PHENO_FILE] \
    --pheno-name [PHENO_COL] \
    --covar [INPUT_COVAR_FILE] \
    --covar-name [COVAR_COL_1],[COVAR_COL_2],...,[COVAR_COL_X] \
    --covar-name sex,EV1,EV6 \
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
    --covar-name [COVAR_COL_1],[COVAR_COL_2],...,[COVAR_COL_X] \
    --meta score \
    --dosage DS \
    --useResidualAsPhenotype \
    --inverseNormal \
    --qtl \
    --out [OUTPUT_PREFIX]
```

