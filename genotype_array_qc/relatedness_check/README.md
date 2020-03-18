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

----

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

----

<details>
<summary>3. Run LD pruning workflow</summary>
</br>

</details>

----

<details>
<summary>4. Identify related individuals to remove before PCA</summary>
</br>

Prior to running PCA, the sample set will be reduced to only unrelated (relatedness degree >3) individuals. This helps mitigate [PC estimation biases caused by high genotypic correlation between samples](https://stats.stackexchange.com/questions/50537/should-one-remove-highly-correlated-variables-before-doing-pca), to produce top PCs that are maximally informative for ancestral diversity. To obtain a subset of unrelated samples, [KING](https://doi.org/10.1093/bioinformatics/btq559) will be used.

<details>
<summary>Standard processing</summary>
</br>

# Standard processing

## Kinship coefficient calculations

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

## Identify unrelated subset of samples

To obtain a list of unrelated individuals, a greedy graph pruning approach is taken.

Sample code:
```R
# TODO: Convert this example code into an R script

library(igraph)
options(stringsAsFactors = F)

# User-specified arguments
input.file <- "king/king.kin0" # KING across-family results table from --kinship
output.file <- "pca_sample_exclusions.txt" # Name for outputting PLINK compatible remove ID file

# Read in KING --kinship output table
king.stats <- read.table(input.file, header = T)

# Combine FID and IIDs into single ID for first component of a pair
edge.heads <- paste0(king.stats$FID1, ":::", king.stats$ID1)
# Combine FID and IIDs into single ID for second component of a pair
edge.tails <- paste0(king.stats$FID2, ":::", king.stats$ID2)
sample.pairs = cbind(edge.heads, edge.tails)

# Create undirected sample graph
sample.graph <- graph.data.frame(sample.pairs, directed = F)

# Data structure to track which samples to exclude
remove.list <- c()

# Get the number of relationships per sample
sample.degrees <- sort(degree(sample.graph), decreasing = T)

# Apply greedy graph pruning approach
# Remove highest degree samples until none are left or
#   the highest degree sample has degree 1
current.graph <- sample.graph
while(length(sample.degrees) > 0 & sample.degrees[1] > 1) {
    remove.list <- c(remove.list, names(sample.degrees[1]))
    current.graph <- current.graph - names(sample.degrees[1])
    sample.degrees = sort(degree(current.graph), decreasing = T)
}

# Update sample pairs post-pruning
unpruned.sample.pairs <- sample.pairs[!(sample.pairs[,1] %in% remove.list) & !(sample.pairs[,2] %in% remove.list),]

# Randomly select one from each pair to remove
random.exclusions <- sapply(1:nrow(unpruned.sample.pairs),
    function(i){unpruned.sample.pairs[i, sample(x = 1:2, size = 1)]})

# Update remove list
remove.list <- c(remove.list, as.vector(random.exclusions))

# Final check that no sample pairs remain
if(sum(!(sample.pairs[,1] %in% remove.list) & !(sample.pairs[,2] %in% remove.list)) == 0){
    stop("Error: Not all sample pairs filtered during graph pruning!")
}

# Export PLINK compatible remove list
final.remove.list <- do.call(rbind, strsplit(remove.list, split = ":::"))
write.table(final.remove.list, file = output.file, 
    sep = " ", row.names = F, col.names = F, quote = F)
```
</details>

<details>
<summary>Parallel processing</summary>
</br>

# Parallel processing

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

----

<details>
<summary>5. Construct unrelated sample set</summary>
</br>

Using the autosome only PLINK file set (prior to LD-pruning), related individuals (degree 3 or less) are removed and variant QC is re-applied.

Sample code:
```shell

# Remove related samples and re-apply variant QC
plink \
    --bfile <bed_bim_fam_prefix> \
    --remove <remove_list_prefix>.txt
    --geno <rate> \
    --hwe <pvalue> \
    --maf <freq> \
    --make-bed \
    --out <output_prefix>
```

Input files:

| FILE | DESCRIPTION |
| --- | --- |
| `<bed_bim_fam_prefix>.bed` | PLINK format bed file for input genotypes |
| `<bed_bim_fam_prefix>.bim` | PLINK format bim file for input genotypes |
| `<bed_bim_fam_prefix>.fam` | PLINK format fam file for input genotypes |
| `<remove_list_prefix>.txt` | Two-column (FID and IID) space-delimited ID file containing samples to remove |

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
| `--bfile <bed_bim_fam_prefix>` | PLINK file set for input genotypes |
| `--remove` | Remove samples based on a list of IDs (FID and IID pairs) |
| `--autosome` | Flag indicating to retain only chromosomes 1-22 |
| `--geno <rate>` | Filters out all variants with missing call rates exceeding `<rate>` (decimal value) |
| `--hwe <pvalue>` | Filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below `<pvalue>` |
| `--maf <freq>` | Filters out all variants with minor allele frequency below `<freq>` (decimal value) |
| `--make-bed` | Flag indicating to generate genotypes in PLINK bed/bim/fam format |
| `--out <output_prefix>` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


----

<details>
<summary>6. Run LD pruning workflow</summary>
</br>

</details>

----

<details>
<summary>7. Run PCA</summary>
</br>

To ensure that calculation of principal components (PCs) will scale reasonably (in terms of both memory and CPU demands) to biobank size data (100k+ samples), [FlashPCA2](https://github.com/gabraham/flashpca) will be used. This PCA algorithm gets its memory and computational efficency from doing block-wise calculations on the data and using the implicitly restarted Arnoldi method in calculating a pre-specificed number of PC approximations. The [manuscript](https://doi.org/10.1093/bioinformatics/btx299) describing FlashPCA2 mechanics highlights that PCA takes under 6 hours using 2GB of RAM for 20 PCs on 500k samples and 100k SNPs.

For the purpose of identifying ancestry-informative SNPs, the PCA will only require calculating the loadings for the top 3 PCs.

</details>

----
