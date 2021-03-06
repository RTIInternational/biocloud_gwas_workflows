# Sex check

## Introduction

This document details the standard analysis workflow for performing a sex check that compares genotypic sex to self-reported sex. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes genotypes in PLINK bed/bim/fam format and a file containing at minimum the FID, IID, and sex for each individual. The workflow produces the following outputs:

1. Output generated by the PLINK `--check-sex` option, including phenotypic sex, genotypic sex and the F estimate for all individuals.
2. File containing only the problematic individuals from file 1.
3. File that can be used with the PLINK `--remove` option to remove individuals that failed the sex check.

## Workflow

The steps in this workflow are as follows:
<details>
<summary>1. Split the X chromosome into PAR and non-PAR</summary>

Sample command:
```
# First merge to account for datasets that have already been split
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --merge-x no-fail \
    --make-bed \
    --out tmp.merge_x

# Now split
plink \
    --bfile tmp.merge_x \
    --split-x b37 \
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
| `--split-x b37 no-fail` | Option telling PLINK to split X based on b37 coordinates and not fail if already split |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>2. LD pruning workflow (separate supporting workflow)</summary>

</details>


<details>
<summary>3. Generate file for providing self-reported sex for sex check</summary>

Sample command:
``` shell
tail -n +[HEADER_ROWS] [PHENO_FILE] |
    perl -lne '
        $delimiter = lc("'[PHENO_DELIMITER]'");	
        $delimiter = ($delimiter eq "comma") ? "," : (($delimiter eq "tab") ? "\t" : (($delimiter eq "space") ? " " : ""));
        chomp;
        @F = split($delimiter)
        print join("\t", $F['[FID_COL]'], $F['[IID_COL]'], $F['[SEX_COL]'])' > [SEX_FILE]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[PHENO_FILE]` | PLINK format bed file for input genotypes |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[SEX_FILE]` | PLINK sex check output for all subjects |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--pheno [PHENO_FILE]` | Phenotype file containing at minimum FID, IID, and sex |
| `--fid_col [FID_COL]` | Column # in phenotype file containing the FID (0-based) |
| `--iid_col [IID_COL]` | Column # in phenotype file containing the IID (0-based) |
| `--sex_col [SEX_COL]` | Column # in phenotype file containing the FID (0-based) |
| `--header_rows [HEADER_ROWS]` | # of header rows in phenotype file |
| `--delimiter [DELIMITER]` | Delimiter used in phenotype file (accepted values comma, tab, space) |
| `--out [SEX_FILE]` | Prefix for output |
</details>


<details>
<summary>(Optional) 4. Generate supplementary freq file from different dataset for sex check for small datasets</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --freqx \
    --out [OUTPUT_PREFIX]
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
| `[OUTPUT_PREFIX].frqx` | Allele frequencies for input genotypes |
| `[OUTPUT_PREFIX].log` | PLINK log file |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bfile [INPUT_BED_BIM_FAM_PREFIX]` | Prefix for input genotypes in PLINK bed/bim/fam format (not the main dataset being QCed) |
| `--out [OUTPUT_PREFIX]` | Prefix for output |
</details>


<details>
<summary>5. Perform sex check</summary>

Sample command:
``` shell
# Small datasets
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --check-sex [FEMALE_MAX_F] [MALE_MIN_F] \
    --update-sex [SEX_FILE] \
    --read-freq [FREQ_FILE] \
    --out [OUTPUT_PREFIX]

# Other datasets
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --check-sex [FEMALE_MAX_F] [MALE_MIN_F] \
    --update-sex [SEX_FILE] \
    --out [OUTPUT_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for input genotypes |
| `[SEX_FILE]` | File generated in step 3 |
| `[FREQ_FILE]` | (Optional) Frequency file for different dataset generated using PLINK `--freq` or `--freqx` command |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUTPUT_PREFIX].sexcheck` | PLINK sex check output for all subjects |
| `[OUTPUT_PREFIX].sexcheck.log` | PLINK log file for sex check |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bfile [INPUT_BED_BIM_FAM_PREFIX]` | Prefix for input genotypes in PLINK bed/bim/fam format |
| `--check-sex [FEMALE_MAX_F] [MALE_MIN_F]` | Flag indicating that PLINK shoud perform a sex check. Default values for `[FEMALE_MAX_F]` and `[MALE_MIN_F]` are 0.2 and 0.8, respectively. |
| `--out [OUTPUT_PREFIX]` | Prefix for output |
</details>


<details>
<summary>6. Generate final files</summary>

Sample command:
``` shell
# Rename output file
perl -lane 'print join("\t",@F);' [SEX_CHECK_RESULTS] > [OUTPUT_PREFIX].sexcheck.all.tsv

# Extract subjects not passing sex check
head -n 1 [OUTPUT_PREFIX].sexcheck.all.tsv > [OUTPUT_PREFIX].sexcheck.problems.tsv
grep PROBLEM [OUTPUT_PREFIX].sexcheck.all.tsv >> [OUTPUT_PREFIX].sexcheck.problems.tsv

# Create remove list
tail -n +2 [OUTPUT_PREFIX].sexcheck.problems.tsv |
    perl -lane 'print join("\t", $F[0], $F[1]);' > [OUTPUT_PREFIX].sexcheck.remove.tsv
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[SEX_CHECK_RESULTS]` | PLINK output from sex check in previous step |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUTPUT_PREFIX].sexcheck.all.tsv` | PLINK sex check output for all subjects |
| `[OUTPUT_PREFIX].sexcheck.problems.tsv` | PLINK sex check output for subjects not passing sex check |
| `[OUTPUT_PREFIX].sexcheck.remove.tsv` | List of subjects not passing sex check that can be fed into PLINK to remove the subjects |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--in [SEX_CHECK_RESULTS]` | Prefix for input genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_PREFIX]` | Prefix for output |
</details>
