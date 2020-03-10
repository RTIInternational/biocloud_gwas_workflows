# Genotype Array QC

## Introduction

This document details the standard analysis workflow for performing QC data from genotyping arrays. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes plus-strand GRCh37 genotypes in PLINK bed/bim/fam format and produces the following outputs:

1. QCed genotypes in PLINK bed/bim/fam format.
2. Summary of variants and subjects removed/flagged during each step of the QC pipeline.

The input and output formats are fully described in the appendix of this document.

## Workflow

The steps in this workflow are as follows:
<details>
<summary>1. Split by chromosome</summary>

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

</details>

<details>
<summary>2. Convert variants to IMPUTE2 ID format</summary>

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
</details>


<details>
<summary>3. Remove duplicate IDs (based on call rate)</summary>

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
<summary>4. Merge chromosomes</summary>

Sample command:
``` shell
for prefix in $([PREFIX_LIST]); do
    if [ [FORMAT] == "bed_bim_fam" ]; then
        echo $prefix.bed $prefix.bim $prefix.fam
    elif [ [FORMAT] == "ped_map" ]; then
        echo $prefix.ped $prefix.map
    fi
done > $fileMergeList

plink \
    --merge-list $fileMergeList \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]

rm $fileMergeList
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |


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
| `--prefix_list [PREFIX_LIST]` | List of prefixes of files to be merged |
| `--format [FORMAT]` | Format of files to be merged (bed_bim_fam, ped_map) |
</details>


<details>
<summary>5. Flag individuals missing chrX or other chromosome</summary>

</details>


<details>
<summary>6. Remove phenotype info in FAM file</summary>

Sample command:
```
perl -pe 's/\S+$/0/;' [INPUT_FAM_FILE]
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
<summary>7. Format phenotype data to standard format</summary>

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
<summary>8. Structure workflow (separate supporting workflow)</summary>

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
<summary>9. Partition data by ancestry</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \ \
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
<summary>10. Call rate filter</summary>

Sample command:
``` shell
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \ \
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
<summary>11. HWE filter</summary>

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
<summary>12. Subject call rate filter (based on autosomes)</summary>

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
<summary>13. Relatedness workflow (separate supporting workflow)</summary>

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
<summary>14. Remove samples based on relatedness</summary>

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
<summary>15. Sex check and sample removal</summary>

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
<summary>16. Excessive homozygosity filtering</summary>

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
<summary>17. Set het haploids to missing</summary>

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


