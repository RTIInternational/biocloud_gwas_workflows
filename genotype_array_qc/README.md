# Genotype Array QC

## Introduction

This document details the standard analysis workflow for performing QC data from genotyping arrays. An automated pipeline, developed using WDL, Cromwell, and Docker, is available for this workflow.

This workflow takes plus-strand GRCh37 genotypes in PLINK bed/bim/fam format. The fam file should include the sex for each individual. The workflow produces the following outputs:

1. QCed genotypes in PLINK bed/bim/fam format.
2. Summary of variants and subjects removed/flagged during each step of the QC pipeline.

The input and output formats are fully described in the appendix of this document.

## Workflow

The steps in this workflow are as follows:
<details>
<summary>1. Split the X chromosome into PAR and non-PAR</summary>

Sample command:
```
plink \
    --bfile [INPUT_BED_BIM_FAM_PREFIX] \
    --split-x b37 no-fail \
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
<summary>2. Remove phenotype info in FAM file</summary>

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
<summary>3. Split by chromosome</summary>

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
<summary>4. Convert variants to IMPUTE2 ID format</summary>

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
<summary>5. Remove duplicate variant IDs (based on call rate) and merge chromosomes</summary>

Sample command:
``` shell
# Get sorted list of all variants
for bim in $(perl -lane 'print $F[1];' [MERGE_LIST]); do
    perl -lane 'print $F[1];' $bim
done | sort > tmp.variants.sorted

# Get sorted list of unique variants
sort -u tmp.variants.sorted > tmp.variants.sorted.unique

# Get list of duplicate variants
comm -23 tmp.variants.sorted \
    tmp.variants.sorted.unique |
    sort -u > \
    tmp.variants.duplicates
    
# Append ___X (where X is a unique number for the variant) to the end of the variant IDs for duplicates
for bim in $(perl -lane 'print $F[1];' [MERGE_LIST]); do
    perl -i.bak -lane '
        BEGIN {
            %duplicates = ();
            open(DUPLICATES, "'tmp'.variants.duplicates");
            while (<DUPLICATES>) {
                chomp;
                $duplicates{$_} = 1;
            }
            close DUPLICATES;
        }
        if (exists($duplicates{$F[1]})) {
            $F[1] = $F[1]."___".($duplicates{$F[1]}++);
        }
        print join("\t", @F);' $bim
done

# Merge chromosomes
/shared/bioinformatics/software/third_party/plink-1.90-beta-4.10-x86_64/plink \
    --merge-list [MERGE_LIST] \
    --make-bed \
    --out tmp.with_duplicates

# Get list of duplicate SNPs
grep ___ tmp.with_duplicates.bim |
    perl -lane 'print $F[1];' > tmp.with_duplicates.duplicates

# Get call rates for duplicate SNPs
/shared/bioinformatics/software/third_party/plink-1.90-beta-4.10-x86_64/plink \
    --bfile tmp.with_duplicates \
    --extract tmp.with_duplicates.duplicates \
    --missing \
    --out tmp.with_duplicates.missing

# Create remove list containing duplicate with higher missing rate
tail -n +2 tmp.with_duplicates.missing.lmiss |
    perl -lane '
      BEGIN {
          %missingness = ();
          %extract = ();
          @variants = ();
      }
      push(@variants, $F[1]);
      $F[1] =~ /^(\S+)___\d+/;
      $baseName = $1;
      if (exists($missingness{$baseName})) {
          if ($F[4] < $missingness{$baseName}) {
              $missingness{$baseName} = $F[4];
              $extract{$baseName} = $F[1];
          }
      } else {
          $missingness{$baseName} = $F[4];
          $extract{$baseName} = $F[1];
      }
      END {
          foreach $variant (@variants) {
              $variant =~ /^(\S+)___\d+/;
              $baseName = $1;
              if ($extract{$baseName} ne $variant) {
                  print $variant;
              }
          }
      }' > tmp.with_duplicates.remove

# Remove duplicates with higher missing rate
/shared/bioinformatics/software/third_party/plink-1.90-beta-4.10-x86_64/plink \
    --bfile tmp.with_duplicates \
    --exclude tmp.with_duplicates.remove \
    --make-bed \
    --out [OUTPUT_BED_BIM_FAM_PREFIX]

# Remove numbers from duplicate variant IDs
perl -i.bak -lne 's/___\d+//; print;' [OUTPUT_BED_BIM_FAM_PREFIX].bim

# Clean up files
rm tmp*
```

Input Files:

| FILE | DESCRIPTION |
| --- | --- |
| `[MERGE_LIST]` | PLINK format merge-list (see https://www.cog-genomics.org/plink/1.9/data#merge_list) |
| `[BED_FILES]` | Array of bed files to merge in same order as `[BIM_FILES]` and `[FAM_FILES]` |
| `[BIM_FILES]` | Array of bim files to merge in same order as `[BED_FILES]` and `[FAM_FILES]` |
| `[FAM_FILES]` | Array of fam files to merge in same order as `[BED_FILES]` and `[BIM_FILES]` |

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
| `--bed_files [BED_FILES]` | Delimited list of bed files to merge in same order as bim and fam files |
| `--bim_files [BIM_FILES]` | Delimited list of bim files to merge in same order as bed and fam files |
| `--fam_files [FAM_FILES]` | Delimited list of fam files to merge in same order as bed and bim files |
| `--out [OUTPUT_BED_BIM_FAM_PREFIX]` | Prefix for output genotypes in PLINK bed/bim/fam format |
</details>


<details>
<summary>6. Remove subjects with >99% missingness</summary>

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
<summary>7. Structure workflow (separate supporting workflow)</summary>

</details>


<details>
<summary>8. Partition data by ancestry</summary>

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
<summary>9. Call rate filter</summary>

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
<summary>10. HWE filter</summary>

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
<summary>11. Set het haploids to missing</summary>

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
<summary>12. Subject call rate filter (based on autosomes)</summary>

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
<summary>13. Excessive homozygosity filtering</summary>

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
<summary>14. Relatedness workflow (separate supporting workflow)</summary>

</details>


<details>
<summary>15. Sex check workflow (separate supporting workflow)</summary>

</details>


<details>
<summary>16. Remove samples based on relatedness (optional)</summary>

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
<summary>17. Remove samples based on discrepant sex (optional)</summary>

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
<summary>18. Merge the PAR and non-PAR regions of the X chromosome</summary>

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
<summary>19. Flag individuals missing chrX or other chromosome</summary>

</details>


