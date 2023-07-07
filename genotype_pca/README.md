# WDL workflow eigenstrat_smartpca: Principle components analysis 

This repository contains a WDL workflow named eigenstrat_smartpca. This workflow performs Eigenstrat PCA analysis on genotype data, including steps for locating high LD regions, removing high LD regions, LD pruning, merging pruned files, extracting LD variants, renaming BIM and FAM files, running SMARTPCA, and creating the final output file.
<br><br>

## Workflow steps

The `eigenstrat_smartpca` workflow consists of the following steps:

<details>
  <summary><b>Step 1:</b> locate_high_ld_regions</summary>
  
   - Description: This step identifies high LD regions in the genotype data.
   - Inputs:
     - `bimfile`: BIM file containing variant information
     - `docker`: Docker image (Ubuntu 18.04)
     - `cpu`: Number of CPUs to allocate
     - `mem`: Amount of memory to allocate
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
     - `mem`: Amount of memory to allocate
    <br>
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
     - `mem`: Amount of memory to allocate
     <br>
</details>

<details>
<summary><b>Step 4:</b>: merge_pruned</summary>
  
   - Description: This step merges the pruned genotype files.
   - Inputs:
     - `pruned_files`: List of pruned genotype files
     - `docker`: Docker image (Ubuntu 18.04)
     - `cpu`: Number of CPUs to allocate
     - `mem`: Amount of memory to allocate
     <br>
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
     - `mem`: Amount of memory to allocate
<br>
</details>


 
<details>
<summary><b>Step 6:</b> rename_bimfam</summary>

   - Description: This step renames the BIM and FAM files.
   - Inputs:
     - `bimfile`: BIM file to rename
     - `famfile`: FAM file to rename
     - `docker`: Docker image (Plink v1.9)
     - `cpu`: Number of CPUs to allocate
     - `mem`: Amount of memory to allocate  
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
     - `mem`: Amount of memory to allocate  
<br>
</details>


<details>
<summary><b>Step 8:</b> create_final_file</summary>

   - Description: This step creates the final output file.
   - Inputs:
     - `study_name`: Name of the study
     - `ancestry`: Ancestry information
     - `evec_file`: SMARTPCA eigenvectors file
     - `famfile`: FAM file containing sample information
     - `docker`: Docker image (Ubuntu 18.04)
     - `cpu`: Number of CPUs to allocate
     - `mem`: Amount of memory to allocate  
<br>
</details>

<br><br>


## Usage

To use this WDL workflow, follow these steps:

1. Clone the repository to your local machine.
1. Modify the `inputs.json` file to provide the required inputs:
   - `bedfile`: Genotype data in BED format (on S3).
   - `bimfile`: Variant information in BIM format (on S3).
   - `famfile`: Sample information in FAM format (on S3).
   - `study_name`: Name of the study.
   - `ancestry`: Ancestry information.
   - Modify any other desired parameters, such as the Docker images and resource allocation (CPU and memory).
1. Execute the WDL workflow. See the [WDL/Cromwell Guide](https://researchtriangleinstitute.sharepoint.com/sites/OmicsGroup/_layouts/15/Doc.aspx?sourcedoc={a2b17bca-8f68-4450-a563-f80609bd497a}&action=edit&wd=target%28Computing%20Infrastructure.one%7Ca745a153-ea3f-4b6e-8f16-9163bfe64932%2FWDL%5C%2FCromwell%20Guide%7C80665feb-2dbf-481d-92d8-cf8c8e7d30dc%2F%29&wdorigin=703) on Microsoft Teams
1. Once the workflow completes, the final output file (`final_file`) will be generated. You can locate this on S3 at: `s3://rti-cromwell-output/cromwell-execution/eigenstrat_smartpca/<job-id>`

<br><br>

### example code

```bash
# clone repo
home=/home/ubuntu
cd $home
git clone --recurse-submodules https://github.com/RTIInternational/biocloud_gwas_workflows

# modify inputs
vim biocloud_gwas_workflows/genotype_pca/inputs.json

# zip dependencies
zip \
    --exclude=*/var/* \
    --exclude=*.git/* \
    --exclude=*/test/* \
    --exclude=*/.idea/* \
    -r imports.zip \
    biocloud_gwas_workflows/

# Open up a connection or tunnel to the Cromwell server using another terminal tab (or more practically, with the screen terminal multiplexer)
ssh -i ~/.ssh/gwas_rsa -L localhost:8000:localhost:8000 ec2-user@54.146.0.138

# Submit job with cURL to Cromwell server (not from within the Cromwell server)
curl -X POST "http://localhost:8000/api/workflows/v1" -H "accept: application/json" \
    -F "workflowSource=@${home}/biocloud_gwas_workflows/genotype_pca/main.wdl" \
    -F "workflowInputs=@${home}/biocloud_gwas_workflows/genotype_pca/inputs.json" \
    -F "workflowDependencies=@${home}/imports.zip" \
    -F "workflowOptions=@${home}/biocloud_gwas_workflows/workflow_options/spot/0216573.000.001_eric_johnson_hiv_omics.json"
# {"id":"6865f67c-a3f9-49aa-8b27-228edc0179a2","status":"Submitted"}

# record job ID
job=033b8637-0dee-429c-87a9-14650e8b9084

# check status of job
curl -X GET "http://localhost:8000/api/workflows/v1/$job/status" # {"status":"Succeeded","id":"033b8637-0dee-429c-87a9-14650e8b9084"}

# download results JSON
curl -X GET "http://localhost:8000/api/workflows/v1/$job/outputs" > outputs.json

# download files from results JSON to local
docker run -it -v $PWD/:/data rtibiocloud/download_wdl_results_from_json:v1_377bef8 \
    --file /data/outputs.json \
    --aws-access-key AKIA12345 \
    --aws-secret-access-key abcde12345
```

<br><br>

## References
For any questions or comments, send me an email or slack and I'll be happy to help: Jesse Marks (jmarks@rti.org)
