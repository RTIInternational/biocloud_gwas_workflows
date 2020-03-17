# Relatedness check

## Introduction

This document details the standard analysis workflow for performing a relatedness check that identifies genetic relatedness using two concepts, [identity-by-state (IBS)](https://isogg.org/wiki/Identical_by_state) and [identity-by-descent (IBD)](https://isogg.org/wiki/Identical_by_descent). An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

As input, this workflow takes genotypes in PLINK bed/bim/fam format. The workflow produces the following outputs:

* A list of IDs corresponding to unrelated samples to keep.
* A list of IDs to remove based on identified duplicates/monozygotic twins.
* A list of IDs to remove based on the [maximum independent vertex sets](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) of identified relatives.
* Pairwise kinship coefficient estimates for all duplicate and related samples up to a user-specific degree.

## Workflow

The steps in this workflow are as follows:

<details>
<summary>1. Remove pedigree info from FAM file </summary>
</details>

In this step, the input PLINK [fam file](https://www.cog-genomics.org/plink/1.9/formats#fam) is modified so that the family ID is unique and parent relationships are removed. An ID map is also created so that PLINK `--update-ids` can be used at a later point to revert the IDs.

Sample command:
```shell
# Create new fam file with pedigree information removed
awk '{$1=NR; $3=0; $4=0; print}' <fam_file> \
    > <updated_fam_file_prefix>.fam 

# Create ID mapping between new IDs and original IDs
awk '{print NR,$2,$1,$2}' <fam_file> \
    > <id_mapping_file_prefix>.txt 
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `<fam_file>` | PLINK format fam file for input genotypes |

Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `<updated_fam_file_prefix>.fam` | PLINK format fam file for output genotypes. Pedigree information removed. |
| `<id_mapping_file_prefix>.txt` | PLINK ID mapping file to go from new IDs to original IDs. Compatible with PLINK `--update-ids` |

<details>
<summary>1. Retain only autosomes</summary>

Sample command:
```
# First merge to account for datasets that have already been split
plink \
    --bfile <bed/bim/fam file prefix> \
    --autosome \
    --make-bed \
    --out <output bed/bim/fam file prefix>
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `<bed/bim/fam file prefix>.bed` | PLINK format bed file for input genotypes |
| `<bed/bim/fam file prefix>.bim` | PLINK format bim file for input genotypes |
| `<bed/bim/fam file prefix>.fam` | PLINK format fam file for input genotypes |


Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `<output bed/bim/fam file prefix>.bed` | PLINK format bed file for output genotypes |
| `<output bed/bim/fam file prefix>.bim` | PLINK format bim file for output genotypes |
| `<output bed/bim/fam file prefix>.fam` | PLINK format fam file for output genotypes |
| `<output bed/bim/fam file prefix>.log` | PLINK log file |


Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bfile <bed/bim/fam file prefix>` | Prefix for input genotypes in PLINK bed/bim/fam format |
| `--autosome` | Flag indicating to retain only chromosomes 1-22 |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out <output bed/bim/fam file prefix>` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>

