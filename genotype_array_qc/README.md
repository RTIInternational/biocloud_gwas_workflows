# Genotype Array QC

## Introduction

This document details the standard analysis workflow for performing QC data from genotyping arrays. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow performs standard QC on genotypes in PLINK bed/bim/fam format. A full description of the inputs and outputs is provided in the appendix of this document.

## Workflow

The steps in this workflow are as follows:
<details>
<summary>Convert to GRCh37 (separate workflow)</summary>

</details>


<details>
<summary>Flip to plus strand (separate workflow)</summary>

</details>



<details>
<summary>Remove phenotype info in FAM file</summary>

Sample command:
```
perl -lane 'print join("\t", @F[0 .. 3], "0\t0");' [INPUT_FAM_FILE]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_FAM_FILE]` | Input FAM file to remove phenotype info from |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[OUTPUT_FAM_FILE]` | Output FAM file phenotype info removed |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--in_fam [INPUT_FAM_FILE]` | Input FAM file to remove phenotype info from |
| `--out_fam [OUTPUT_FAM_FILE]` | Output FAM file phenotype info removed |
</details>


<details>
<summary>Convert variants to IMPUTE2 ID format (separate workflow)</summary>

</details>


<details>
<summary>Remove subjects with >99% missingness</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --mind 0.99 \
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
| `--mind 0.99` | Option indicating that individuals with >99% missingness should be excluded |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Structure workflow (separate workflow)</summary>

</details>


<details>
<summary>Partition data by ancestry</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --keep [KEEP_LIST] \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for input genotypes |
| `[KEEP_LIST]` | List of subjects to keep |


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
| `--keep [KEEP_LIST]` | List of subjects to keep |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Call rate filter</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --geno [CALL_RATE_THRESHOLD] \
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
| `--geno [CALL_RATE_THRESHOLD]` | Call rate threshold for excluding SNPs (e.g., 0.01) |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>HWE filter</summary>

Sample command:
``` shell
# Autosomes
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --hwe [HW_PVALUE_THRESHOLD] \
    --autosome \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX].autosome_hwe

# Chr X for females
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --hwe [HW_PVALUE_THRESHOLD] \
    --filter-females \
    --chr 23 \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe_females

# Extract chr X SNPs from full dataset
perl -lane 'print $F[1];' [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe_females.bim > \
    [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe_females.extract
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --extract [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe_females.extract \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe

# Merge autosomes and chr X
plink \
  --bfile [OUTPUT_BED_BIM_FAM_PREFIX].autosome_hwe \
  --bmerge [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe.bed \
          [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe.bim \
          [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe.fam \
  --make-bed \
  --out [OUTPUT_BED_BIM_FAM_PREFIX]

# Remove intermediates
rm [OUTPUT_BED_BIM_FAM_PREFIX].autosome_hwe*
rm [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe_females*
rm [OUTPUT_BED_BIM_FAM_PREFIX].x_hwe*
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
| `--hwe [HW_PVALUE_THRESHOLD]` | HW p-value threshold for excluding SNPs (e.g., 0.0001) |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Set het haploids to missing</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --set-hh-missing \
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
| `--set-hh-missing` | Flag indicating that PLINK should set heterozygous haploids to missing |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Subject call rate filter (based on autosomes)</summary>

Sample command:
``` shell
# Autosomes
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --mind [CALL_RATE_THRESHOLD] \
    --autosome \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX].autosomes

# Extract subjects from full dataset
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --keep [OUTPUT_BED_BIM_FAM_PREFIX].autosomes.fam \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]

# Remove intermediates
rm [OUTPUT_BED_BIM_FAM_PREFIX].autosomes*
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
| `--mind [CALL_RATE_THRESHOLD]` | Call rate threshold for excluding subjects (e.g., 0.01) |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Excessive homozygosity filtering</summary>

Sample command:
``` shell
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
</details>


<details>
<summary>Relatedness workflow (separate workflow)</summary>

</details>


<details>
<summary>Sex check workflow (separate workflow)</summary>

</details>


<details>
<summary>Remove samples based on relatedness (optional)</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --remove [REMOVE_LIST] \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for input genotypes |
| `[REMOVE_LIST]` | List of subjects to remove |


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
| `--remove [REMOVE_LIST]` | List of subjects to remove |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Remove samples based on discrepant sex (optional)</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --remove [REMOVE_LIST] \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[INPUT_BED_BIM_FAM_PREFIX].bed` | PLINK format bed file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].bim` | PLINK format bim file for input genotypes |
| `[INPUT_BED_BIM_FAM_PREFIX].fam` | PLINK format fam file for input genotypes |
| `[REMOVE_LIST]` | List of subjects to remove |


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
| `--remove [REMOVE_LIST]` | List of subjects to remove |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>Merge the PAR and non-PAR regions of the X chromosome</summary>

Sample command:
```
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --merge-x no-fail \
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
| `--merge-x no-fail` | Option telling PLINK to merge the PAR and non-PAR regions and not fail if already split |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>



<details>
<summary>Flag individuals missing chrX or other chromosome</summary>

</details>


