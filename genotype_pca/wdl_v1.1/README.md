# WDL workflow eigenstrat_smartpca: Principle components analysis 

This repository contains a WDL workflow named eigenstrat_smartpca. This workflow performs Eigenstrat PCA analysis on Plink formatted genotype data (bed,bim,fam). It uses a user-provided [UCSC Genome Browser's BED formatted](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file to locate regions of high LD and subsequently remove SNP from the Plink files that are within those genomic loci. The workflow then runs smartpca – a program implemented within the software package eigenstrat that performs PCA on genetic data to identify population structure and correct for population stratification in association studies. It calculates principal components and eigenvalues from the genetic data matrix. The final result is a file that contains the principal component scores for each individual (samples) in your dataset along the corresponding eigenvectors. Each row represents one individual, and the columns contain the scores for each principal component. 

### example final evec results file:
```
           #eigvals:     1.595     1.211     1.148     1.105     1.098     1.096     1.092     1.086     1.079     1.078
           ID_1:ID_1     0.0502      0.0046      0.0011     -0.0209      0.0739      0.0050     -0.0314      0.0452      0.0056      0.0564              ???
           ID_2:ID_2    -0.1909     -0.0023     -0.1874      0.1941     -0.0015     -0.2971     -0.4311      0.3234     -0.2767     -0.0264              ???
           ID_3:ID_3     0.0472     -0.0227      0.0344      0.0158      0.0155      0.0322     -0.0678      0.0004     -0.0423     -0.0128              ???
```
**_Notes_**
* The header of the final results file includes the eigenvalues (eigvals).
* We do not include case/control information, thus the last column will contain `???`. Otherwise it would have indicated case or control.
* SampleIDs were renamed to avoid program crashes caused by length issues, but their order matches the provided Plink genotype data.
<br>


## EIGENSOFT
[EIGENSOFT](https://www.hsph.harvard.edu/alkes-price/software/) is a widely used software package for population genetics and principal component analysis (PCA). It was developed by the Population Genetics Working Group at the Department of Genetics, Harvard Medical School. Eigensoft provides various tools for analyzing genetic data, particularly single-nucleotide polymorphism (SNP) data obtained from genotyping arrays or sequencing studies. 

Eigensoft's primary functionality is conducting PCA on large-scale genetic data. PCA is used to visualize and identify population structure and genetic relatedness patterns within a dataset. It helps researchers understand the genetic diversity and relationships among individuals or populations. **Eigenstrat** is a part of the Eigensoft software package that employs PCA (in particular smartpca). 

### Detecting Population Stratification
Eigenstrat employs principal component analysis (PCA) on genetic data to identify and correct for population stratification. Here's how it typically works:
<details>
  <summary>expand for explanation</summary>

1. **Genetic Data Input**: Eigenstrat takes genetic data, usually in the form of SNP (single-nucleotide polymorphism) data, as input. This data contains genotypes of individuals at various genetic markers across the genome.
1. **PCA Calculation**: Eigenstrat performs principal component analysis (PCA) on the genetic data. PCA helps identify the principal components that explain the largest sources of genetic variation within the dataset. In the context of population stratification, these principal components often correspond to the underlying genetic ancestry of the individuals.
1. **Population Structure Inference**: The principal components identified by PCA can be used to infer the genetic ancestry of each individual in the dataset. This information allows Eigenstrat to group individuals with similar genetic backgrounds together.
1. **Statistical Correction**: After inferring the population structure, Eigenstrat statistically corrects for the confounding effects of population stratification in association tests. It adjusts the association statistics to account for the genetic ancestry of the individuals, thereby reducing the risk of false-positive associations.
   
By applying Eigenstrat's population stratification correction, researchers can improve the accuracy and reliability of genetic association studies, making it a valuable tool in genetic epidemiology and related research areas.
</details>

<br>

## Workflow steps
The workflow consists of multiple tasks executed within the `eigenstrat_smartpca_wf.wdl` file. The tasks in this workflow involve executing bash, python, and plink commands to perform various operations. 

<details>
<summary><b>Step 1:</b> locate_high_ld_regions</summary>

  - Description: This step identifies high LD regions in the genotype data.
  - Inputs:
    - `bimfile`: BIM file containing variant information
    - `reference_file`: `Tab Separated text file containing regions of high LD.` See https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD) for examples.
    - `docker`: Docker image (Ubuntu 22.04)
    - `cpu`: Number of CPUs to allocate
    - `mem_gb`: Amount of memory to allocate
</details>

<details>
<summary><b>Step 2:</b> remove_high_ld_regions</summary>
  
   - Description: This step removes high LD regions from the genotype data.
   - Inputs:
     - `study_name`: Name of the study
     - `ancestry`: Ancestry information
     - `bedfile`: BED file containing genotype data
     - `bimfile`: BIM file containing variant information
     - `famfile`: FAM file containing sample information
     - `high_ld_regions`: High LD regions identified in the previous step
     - `docker`: Docker image (Plink v1.9)
     - `cpu`: Number of CPUs to allocate
     - `mem_gb`: Amount of memory to allocate
</details>

<details>
<summary><b>Step 3:</b> ld_pruning</summary>

   - Description: This step performs LD pruning on the genotype data.
   - Inputs:
     - `study_name`: Name of the study
     - `ancestry`: Ancestry information
     - `bedfile`: BED file containing genotype data
     - `bimfile`: BIM file containing variant information
     - `famfile`: FAM file containing sample information
     - `docker`: Docker image (Plink v1.9)
     - `cpu`: Number of CPUs to allocate
     - `mem_gb`: Amount of memory to allocate
</details>

<details>
<summary><b>Step 4:</b>: merge_pruned</summary>
  
   - Description: This step merges the pruned genotype files.
   - Inputs:
     - `pruned_files`: List of pruned genotype files
     - `docker`: Docker image (Ubuntu 18.04)
     - `cpu`: Number of CPUs to allocate
     - `mem_gb`: Amount of memory to allocate
</details>

<details>
<summary><b>Step 5:</b> extract_ld_variants</summary>

  - Description: This step extracts LD variants from the genotype data.
   - Inputs:
     - `study_name`: Name of the study
     - `ancestry`: Ancestry information
     - `bedfile`: BED file containing genotype data
     - `bimfile`: BIM file containing variant information
     - `famfile`: FAM file containing sample information
     - `combined_variants`: Combined variant information from the previous step
     - `docker`: Docker image (Plink v1.9)
     - `cpu`: Number of CPUs to allocate
     - `mem_gb`: Amount of memory to allocate
</details>

<details>
<summary><b>Step 6:</b> rename_bimfam</summary>

   - Description: This step renames the BIM and FAM files.
   - Inputs:
     - `bimfile`: BIM file to rename
     - `famfile`: FAM file to rename
     - `docker`: Docker image (Plink v1.9)
     - `cpu`: Number of CPUs to allocate
     - `mem_gb`: Amount of memory to allocate  
</details>

<details>
<summary><b>Step 7:</b> run_smartpca</summary>

   - Description: This step runs SMARTPCA analysis on the genotype data.
   - Inputs:
     - `ancestry`: Ancestry information
     - `study_name`: Name of the study
     - `bedfile`: BED file containing genotype data
     - `bimfile`: Renamed BIM file
     - `famfile`: Renamed FAM file
     - `docker`: Docker image (Eigensoft v6.1.4)
     - `cpu`: Number of CPUs to allocate
     - `mem_gb`: Amount of memory to allocate  
</details>
<br>

## Creation of AWS Healthomics Workflow

``` bash
aws_access_key_id="<AWS_ACCESS_KEY_ID>"
aws_secret_access_key="<AWS_SECRET_ACCESS_KEY>"
aws_region_name="<AWS_REGION_NAME>"
repo_dir="<REPO_DIR>"
docker run -ti \
    -v $repo_dir:$repo_dir \
    -e task=create_wf \
    -e aws_access_key_id=$aws_access_key_id \
    -e aws_secret_access_key=$aws_secret_access_key \
    -e aws_region_name=$aws_region_name \
    -e repo_dir=$repo_dir \
    -e main=$repo_dir/genotype_pca/wdl_v1.1/eigenstrat_smartpca_wf.wdl \
    -e parameter_template=$repo_dir/genotype_pca/wdl_v1.1/eigenstrat_smartpca_parameters.json \
    -e name="eigenstrat_smartpca" \
    -e description="Workflow for calculating genotype PCs from plink binary format files" \
    -e engine=WDL \
    -e storage_capacity=2000 \
    --rm rtibiocloud/healthomics_tools:v1.0_f456174
```
Note: creation of workflow in AWS Healthomics only needs to be done before the first use and upon modifications to the workflow.
<br>


## Usage

To use this WDL workflow, follow these steps:

1. Modify the `config_templates/eigenstrat_smartpca_wf_args.json` file to provide the required inputs:
   - `study_name`: Name of the study.
   - `ancestry`: Ancestry information.
   - `bedfile`: Genotype data in BED format (on S3).
   - `bimfile`: Variant information in BIM format (on S3).
   - `famfile`: Sample information in FAM format (on S3).
   - Modify any other desired parameters, such as the Docker images and resource allocation (CPU and memory).
2. Execute the workflow using AWS Healthomics
``` bash
aws_access_key_id="<AWS_ACCESS_KEY_ID>"
aws_secret_access_key="<AWS_SECRET_ACCESS_KEY>"
aws_region_name="<AWS_REGION_NAME>"
docker run -ti \
    -v <HOST_DIR>:<CONTAINER_DIR> \
    -e task=start_run \
    -e charge_code=<CHARGE_CODE> \
    -e aws_access_key_id=<AWS_ACCESS_KEY_ID> \
    -e aws_secret_access_key=<AWS_SECRET_ACCESS_KEY> \
    -e aws_region_name=<AWS_REGION_NAME> \
    -e workflow_id=<WORKFLOW_ID> \
    -e parameters=<PARAMETERS> \
    -e name=<NAME> \
    -e output_uri=<OUTPUT_URI> \
    -e run_metadata_output_dir=<RUN_METADATA_OUTPUT_DIR> \
    -e storage_capacity=<STORAGE_CAPACITY> \
    --rm rtibiocloud/healthomics_tools:v1.0_f456174

```

<br>

## regions of high LD
Provide [BED-formatted](https://en.wikipedia.org/wiki/BED_(file_format)) files. BED starts are zero-based and BED ends are one-based ([ref](https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/bedtools/BEDTools-User-Manual.v4.pdf#%5B%7B%22num%22%3A381%2C%22gen%22%3A0%7D%2C%7B%22name%22%3A%22XYZ%22%7D%2C54%2C392.724%2C0%5D)). The workflow handles this with
```
if region[1] < pos <= region[2]:
```
<details>
  <summary>example BED files</summary>

There are example BED files on S3 for genome build 37 and 38. These files were created using the wiki at: https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD).
You can use these, or create your own to provide to the workflow.
- `s3://rti-bioinformatics-resources/linkage_disequilibrium/regions_of_high_ld_for_pca_wdl_wf_hg19.bed`
- `s3://rti-bioinformatics-resources/linkage_disequilibrium/regions_of_high_ld_for_pca_wdl_wf_hg38.bed`
</details>

<br>

## Contact
For any questions or comments, send me an email and I'll be happy to help: Jesse Marks (jmarks@rti.org)
