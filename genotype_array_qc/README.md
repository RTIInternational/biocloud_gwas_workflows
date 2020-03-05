# Genotype Array QC

## Introduction

This document details the standard analysis workflow for performing QC data from genotyping arrays. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes plus-strand GRCh37 genotypes in PLINK bed/bim/fam format and produces the following outputs:

1. QCed genotypes in PLINK bed/bim/fam format.
2. Summary of variants and subjects removed/flagged during each step of the QC pipeline.

The input and output formats are fully described in the appendix of this document.

The steps in this workflow are as follows:

1. Split by chromosome
2. Convert variants to IMPUTE2 ID format
3. Remove duplicate IDs (based on call rate)
4. Merge chromosomes
5. Flag individuals missing chrX or other chromosome
6. Remove phenotype info in FAM file
7. Format phenotype data to standard format
8. Structure workflow (separate supporting workflow)
9. Partition data by ancestry
10. Call rate filter
11. HWE filter
12. Subject call rate filter (based on autosomes)
13. Relatedness workflow (separate supporting workflow)
14. Remove samples based on relatedness
15. Sex check and sample removal
16. Excessive homozygosity filtering
17. Set het haploids to missing

Each of these steps in described in detail below.

### 1. Split by chromosome

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --chr [CHR] \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]
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
| `[OUTPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for output genotypes |
| `[OUTPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for output genotypes |
| `[OUTPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for output genotypes |
| `[OUTPUT_BED_BIM_FAM_PREFIX].log` | PLINK log file |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bfile [INPUT_BED_BIM_FAM_PREFIX]` | Prefix for input genotypes in PLINK bed/bim/fam format |
| `--chr [CHR]` | Chromosome to extract (1-26, X, Y, XY, MT) |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |


### 2. Convert variants to IMPUTE2 ID format

Sample command:
``` shell
convert_to_1000g_ids.pl \
    --file_in [INPUT_BIM_FILE] \
    --file_out [OUTPUT_BIM_FILE] \
    --legend [INPUT_1000G_LEGEND_FILE] \
    --file_in_id_col [ID_COL_NUM] \
    --file_in_chr_col [CHR_COL_NUM] \
    --file_in_pos_col [POS_COL_NUM] \
    --file_in_a1_col [A1_COL_NUM] \
    --file_in_a2_col [A2_COL_NUM] \
    --chr [CHR]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BIM_FILE]` | PLINK format bim file |
| `[INPUT_1000G_LEGEND_FILE]` | IMPUTE2 1000G legend file |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUTPUT_BIM_FILE]` | PLINK format bim file with IDs in IMPUTE2 format |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--file_in [INPUT_BIM_FILE]` | Path of input bim file |
| `--file_out [OUTPUT_BIM_FILE]` | Path of output bim file |
| `--legend [INPUT_1000G_LEGEND_FILE]` | Path of IMPUTE2 1000G legend file |
| `--file_in_id_col [ID_COL_NUM]` | ID column number (zero-based) |
| `--file_in_chr_col [CHR_COL_NUM]` | Chromosome column number (zero-based) |
| `--file_in_pos_col [POS_COL_NUM]` | Position column number (zero-based) |
| `--file_in_a1_col [A1_COL_NUM]` | Allele 1 column number (zero-based) |
| `--file_in_a2_col [A2_COL_NUM]` | Allele 2 column number (zero-based) |
| `--chr [CHR]` | Chromosome (1-22, X_NONPAR, PAR1, PAR2) |
