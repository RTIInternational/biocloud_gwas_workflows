# LD Pruning

## Introduction

This document details the standard analysis workflow for performing LD pruning of genotypes. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes the following inputs:
1. Genotypes in PLINK bed/bim/fam format
2. Values for window size, step size, and r^2 threshold for pruning (see https://www.cog-genomics.org/plink/1.9/ld)

This workflow generates the following outputs:
1. LD-pruned genotypes in bed/bim/fam format

## Workflow

The steps in this workflow are as follows:
<details>
<summary>1. Generate pruning lists</summary>

Sample command:
```
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --indep-pairwise [WINDOW_SIZE] [STEP_SIZE] [RSQ_THRESHOLD] \
    --out [OUT_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for input genotypes |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUT_PREFIX].prune.in` | List of variants in approximate linkage equilibrium |
| `[OUT_PREFIX].prune.out` | List of excluded variants |
| `[OUT_PREFIX].log` | PLINK log file |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bfile [INPUT_BED_BIM_FAM_PREFIX]` | Prefix for input genotypes in PLINK bed/bim/fam format |
| `--window_size [WINDOW_SIZE]` | Window size to use for call to `--indep-pairwise` |
| `--step_size [STEP_SIZE]` | Step size to use for call to `--indep-pairwise` |
| `--rsq_threshold [RSQ_THRESHOLD]` | R-squared threshold to use for call to `--indep-pairwise` |
| `--out [OUTPUT_PREFIX]` | Prefix for output files |
</details>


<details>
<summary>2. Generate LD pruned genotype files</summary>

Sample command:
```
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --extract [PRUNE_IN_FILE] \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for input genotypes |
| `[PRUNE_IN_FILE]` | List of variants in approximate linkage equilibrium from step 1 |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUTPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for output genotypes |
| `[OUTPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for output genotypes |
| `[OUTPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for output genotypes |
| `[OUTPUT_BED_BIM_FAM_PREFIX].log` | PLINK log file |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bfile [INPUT_BED_BIM_FAM_PREFIX]` | Prefix for input genotypes in PLINK bed/bim/fam format |
| `--prune_in [PRUNE_IN_FILE]` | List of variants in approximate linkage equilibrium from step 1 |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>
