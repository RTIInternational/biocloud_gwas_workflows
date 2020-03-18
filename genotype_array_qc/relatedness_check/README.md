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
</br>

In this step, the input PLINK [fam file](https://www.cog-genomics.org/plink/1.9/formats#fam) is modified so that the family ID is unique and parent relationships are removed. An ID map is also created so that PLINK `--update-ids` can be used at a later point to revert the IDs.

Sample command:
```shell
# Create new fam file with pedigree information removed
# Column 1 is uniquely assigned as the corresponding row number
# Columns 3 and 4 are assigned all 0s
awk '{$1=NR; $3=0; $4=0; print}' <fam_file> \
    > <updated_fam_file_prefix>.fam 

# Create ID mapping between new IDs and original IDs
# Column order = new FID, new IID, old FID, old IID
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

</details>

<details>
<summary>2. Filter to only autosomes and perform initial quality filtering</summary>
</br>

Only the autosomes are needed for accurate calculations of relatedness, so all other chromosomes are excluded to simplify the workflow and reduce resource requirements. For each variant, the genotype call rate, minor allele frequency, and Hardy-Weinberg equilibrium are used to filter out low-quality variants that would bias downstream calculations. 

Sample command:
```shell
# Apply variant filtering and reduce data set to only autosomes
# Explicit input specification of bim/bed/fam separately allows 
#   for incorporation of updated fam file without pedigree info
plink \
    --bed <bed_prefix>.bed \
    --bim <bim_prefix>.bim \
    --fam <fam_prefix>.fam \
    --autosome \
    --geno <rate> \
    --hwe <pvalue> \
    --maf <freq> \
    --make-bed \
    --out <output_prefix>
```

Input files:

| FILE | DESCRIPTION |
| --- | --- |
| `<bed_prefix>.bed` | PLINK format bed file for input genotypes |
| `<bim_prefix>.bim` | PLINK format bim file for input genotypes |
| `<fam_prefix>.fam` | PLINK format fam file for input genotypes. Pedigree information removed. |

Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `<output_prefix>.bed` | PLINK format bed file for output genotypes |
| `<output_prefix>.bim` | PLINK format bim file for output genotypes |
| `<output_prefix>.fam` | PLINK format fam file for output genotypes |
| `<output_prefix>.log` | PLINK log file |

Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--bed <bed_prefix>.bed` | PLINK bed file for input genotypes |
| `--bim <bim_prefix>.bed` | PLINK bim file for input genotypes |
| `--fam <fam_prefix>.bed` | PLINK fam file for input genotypes |
| `--autosome` | Flag indicating to retain only chromosomes 1-22 |
| `--geno <rate>` | Filters out all variants with missing call rates exceeding `<rate>` (decimal value) |
| `--hwe <pvalue>` | Filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below `<pvalue>` |
| `--maf <freq>` | Filters out all variants with minor allele frequency below `<freq>` (decimal value) |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out <output_prefix>` | Prefix for output genotypes in PLINK bed/bim/fam format |

</details>

<details>
<summary>3. Run LD pruning workflow</summary>
</br>

</details>

<details>
<summary>4. Exclude related individuals for PCA</summary>
</br>

Prior to running PCA, the sample set will be reduced to only unrelated (relatedness degree >3) individuals. This helps mitigate [PC estimation biases caused by high genotypic correlation between samples](https://stats.stackexchange.com/questions/50537/should-one-remove-highly-correlated-variables-before-doing-pca), to produce top PCs that are maximally informative for ancestral diversity. To obtain a subset of unrelated samples, [KING](https://doi.org/10.1093/bioinformatics/btq559) will be used.

<details>
<summary>Standard processing</summary>
</br>

Sample command:
```shell
# Calculate kinship coefficients and report only pairwise 
#   relationships greater than or equal to the specified degree
king \
    --cpus <num_cpus> \
    -b <bed_prefix>.bed \
    --fam <fam_prefix>.fam \
    --bim <bim_prefix>.bim \
    --kinship \
    --degree <degree> \
    --prefix <output_prefix>
```

Input files:

| FILE | DESCRIPTION |
| --- | --- |
| `<bed_prefix>.bed` | PLINK format bed file for input genotypes |
| `<bim_prefix>.bim` | PLINK format bim file for input genotypes |
| `<fam_prefix>.fam` | PLINK format fam file for input genotypes |

Output Files:

| FILE | DESCRIPTION |
| --- | --- |
| `<output_prefix>.kin` | Within-family kinship calculation results table. For PLINK fam files with no pedigree information, this file should only contain a header line with no table records. |
| `<output_prefix>.kin0` | Across family kinship calculation results table. Column descriptions [here](http://people.virginia.edu/~wc9c/KING/manual.html#WITHIN) |

Parameters:

| PARAMETER | DESCRIPTION |
| --- | --- |
| `--b <bed_prefix>.bed` | PLINK bed file for input genotypes |
| `--bim <bim_prefix>.bed` | PLINK bim file for input genotypes |
| `--fam <fam_prefix>.bed` | PLINK fam file for input genotypes |
| `--cpus <num_cpus>` | Flag indicating to retain only chromosomes 1-22 |
| `--kinship` | Filters out all variants with minor allele frequency below `<freq>` (decimal value) |
| `--degree <degree>` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--prefix <output_prefix>` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>

<details>
<summary>Parallel processing</summary>
</br>

As descrbined in the KING [documentation](http://people.virginia.edu/~wc9c/KING/manual.html#WITHIN), the `--proj` option can be used to take advantage of batch sample processing which lends itself to easy parallelization. This approach offers both real time computational speedups and reduced memory requirements per run when compared to processing the whole data set at once. Batch sample processing becomes necessary when analyzing tens of thousands of samples or more. For illustrative purposes, assume that a data set is partitioned into 3 batches of samples. Six iterations of KING would need to be run to perform all pairwise comparisons between all samples across the batches.

**Note:** If all the partitions have *exactly* the same sample size, then the value provided to the `--proj` option is simply the sample size of any partition. If the partitions differ in size, then the value provided to the `--proj` option should be the sample size of the *first* partition listed for the `-b`, `--fam`, and `--bim` options.

Sample command:
```shell
# Calculate kinship coefficients between partitions 1 and 2
king \
    --cpus <num_cpus> \
    -b <partition1>.bed,<partition2>.bed \
    --fam <partition1>.fam,<partition2>.fam \
    --bim <partition1>.bim,<partition2>.bim \
    --kinship \
    --degree <degree> \
    --proj <partition1_size> \
    --prefix <partition12_prefix>
    
# Calculate kinship coefficients between partitions 2 and 3
king \
    --cpus <num_cpus> \
    -b <partition2>.bed,<partition3>.bed \
    --fam <partition2>.fam,<partition3>.fam \
    --bim <partition2>.bim,<partition3>.bim \
    --kinship \
    --degree <degree> \
    --proj <partition2_size> \
    --prefix <partition23_prefix>

# Calculate kinship coefficients between partitions 1 and 3
king \
    --cpus <num_cpus> \
    -b <partition1>.bed,<partition3>.bed \
    --fam <partition1>.fam,<partition3>.fam \
    --bim <partition1>.bim,<partition3>.bim \
    --kinship \
    --degree <degree> \
    --proj <partition1_size> \
    --prefix <partition13_prefix>

# Calculate kinship coefficients within each partition
king \
    --cpus <num_cpus> \
    -b <partition1>.bed \
    --fam <partition1>.fam \
    --bim <partition1>.bim \
    --kinship \
    --degree <degree> \
    --prefix <partition1_prefix>

king \
    --cpus <num_cpus> \
    -b <partition2>.bed \
    --fam <partition2>.fam \
    --bim <partition2>.bim \
    --kinship \
    --degree <degree> \
    --prefix <partition2_prefix>

king \
    --cpus <num_cpus> \
    -b <partition3>.bed \
    --fam <partition3>.fam \
    --bim <partition3>.bim \
    --kinship \
    --degree <degree> \
    --prefix <partition3_prefix>

```
</details>

</details>
