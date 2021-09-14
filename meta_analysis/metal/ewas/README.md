# EWAS meta-analysis workflow: METAL
Documentation for running an EWAS meta-analysis using the [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) software tool.

## Pre-requisites
* Unix-based operating system (Linux or OSx. Sorry Windows folks.)

<br><br><br>


## Workflow overview
This workflow performs an EWAS meta-analysis.

INPUT:

A list, of any size, of EWAS summary statistics. THERE MUST BE ONLY ONE SUMMARY STATS FILE PER STUDY. So if your results are split up by chromosome, you must combine those files.

OUTPUT: 
* Meta-analysis results with singletons removed (probe IDs only present in one study). Split by chromosome.
* P-value filtered meta-analysis results.
* Manhattan and QQ plots of meta-analysis results (no singletons).

<br><br><br>

## Running the analysis on AWS Batch
See [this cheat sheet](https://teams.microsoft.com/l/channel/19%3Af42632e48b7c4b9e9f362afa1e4e1957%40thread.tacv2/tab%3A%3A61aecad5-13fa-4bde-adce-ba3b16950439?groupId=9179c917-4161-4094-bec2-b13d4862274c&tenantId=2ffc2ede-4d44-4994-8082-487341fa43fb) on Microsoft Teams.

1. Clone the parent repository biocloud_gwas_workflows to your local machine.
```
cd /home/ubuntu/
git clone https://github.com/RTIInternational/biocloud_gwas_workflows
```

2. Navigate to ewas folder.
```
cd biocloud_gwas_workflows/meta_analysis/metal/ewas/
```

3. Edit the input configuration file located in the `config_templates` folder.

<details>
   <summary>Expand for a line-by-line description of the input.json file</summary>
   
   1. Give your analysis a name (String).
   ```
  "metal_ewas_meta_analysis_wf.plot_basename": "<enter_the_plot_name_here_please>",
   
   ```

   2. Set the p-value threshold to filter (Float).
   ```
   "metal_ewas_meta_analysis_wf.pvalue_threshold": 0.001,
   ```
 
   3. Specify the ancestry group (String).
   ```
      "metal_ewas_meta_analysis_wf.ancestry": "afr",   
   ```
   
   4. Specify which chromosomes to perform the meta-analysis on (Array[Int]).
   ```
     "metal_ewas_meta_analysis_wf.chromosomes_to_keep": [1,2],   
   ```
   
   5. Specify which column the probe ID is in for each study/cohort (Array[Int]).
   ```
     "metal_ewas_meta_analysis_wf.probe_id_column": [1,1],   
   ```
   
   6. Specify which column the chromosome is in for each study/cohort (Array[Int]).
   ```
     "metal_ewas_meta_analysis_wf.chromosome_column": [6,6],   
   ```
   
   7. Specify which column the position is in for each study/cohort (Array[Int]) 
   ```
     "metal_ewas_meta_analysis_wf.position_column": [7,7],   
   ```
   
   8. Specify which column the effect is in for each study/cohort (Array[Int]).
   ```
     "metal_ewas_meta_analysis_wf.effect_size_column": [2,2],   
   ```

   9. Specify which column the standard error is in for each study/cohort (Array[Int]).
   ```
     "metal_ewas_meta_analysis_wf.standard_error_column": [3,3],   
   ```
   
   10. Specify which column the pvalue is in for each study/cohort (Array[Int]).
   ```
     "metal_ewas_meta_analysis_wf.pvalue_column": [4,4],   
   ```
   
   11. Specify a moniker for each study/cohort (Array[String]).
   ```
  "metal_ewas_meta_analysis_wf.study_basename": [
  "study_number_one",
  "study_number_two"
  ],
   ```
   
   12. Specify the AWS S3 locations of the summary statistics for each study/cohort (Array[String]).

  "metal_ewas_meta_analysis_wf.ewas_results": [
  "s3://main-bucket/the/location/of/study_number_one/summary_stats/study_number_one.csv.gz",
  "s3://main-bucket/the/location/of/study_number_two/summary_stats/study_number_two.csv",
  ],
   
   13. Specify if the summary statistics are comma separated for each study/cohort (Array[true or false])
  "metal_ewas_meta_analysis_wf.comma_separated": [
      true,
      false
  ]
 
   </details><br>
   
   4. Make a zip file of the git repo you cloned so Cromwell can handle the local WDL imports.
   ```
   # Change to directory immediately above the biocloud_gwas_workflows repository
   cd /path/to/biocloud_gwas_workflows/
   cd ..
   # Make zipped copy of repo somewhere
   zip --exclude=*var/* --exclude=*.git/* -r /path/to/biocloud_gwas_workflows/meta_analysis/metal/ewas/biocloud_gwas_workflows.zip biocloud_gwas_workflows/
   ```
  
   5. Submit via an API call to a Cromwell server on AWS. Change the paths for each parameter to reflect the correct file paths on your machine. Note to specify the appropriate batch environment so that the correct project gets charged. See the list of charge-code files [here](https://github.com/RTIInternational/bioinformatics/tree/master/config/aws_batch_queues). Contact Jesse Marks (jmarks@rti.org) if you need a new to create a new batch environment for your project.
   ```
   curl -X POST "http://localhost:8000/api/workflows/v1" -H "accept: application/json" \
       -F "workflowSource=@/path/to/biocloud_gwas_workflows/meta_analysis/metal/ewas/full_ewas_meta.wdl" \
       -F "workflowInputs=@/path/to/biocloud_gwas_workflows/meta_analysis/metal/ewas/config_templates/inputs.json" \
       -F "workflowDependencies=@/path/to/biocloud_gwas_workflows/meta_analysis/metal/ewas/biocloud_gwas_workflows.zip"
       -F "workflowOptions=@/path/to/the/charge_code/json_file/project_with_money_charge_code.json"
   
   ### Record the job ID that prints out upon execution of the above command
   ### It will be a string of digits with a format such as `121cecfd-046c-417f-95d3-90225dc833c6`
   ```
   
   6. When the analysis is complete, you can retrieve the results from S3. They will be located in the directory
   ```
   s3://rti-cromwell-output/cromwell-execution/<job-id>/
   ```
   
   <br><br><br>
   
   
   
   
   ## Running the analysis locally
   Sometimes AWS Batch and Cromwell don't play nice together. In that case you can also run the workflow locally. Here are those steps.
   
   ### Prerequisites
   
   * Java v1.8 or higher
   * [Docker](https://docs.docker.com/get-docker/)
   * [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html)
  
1. Install [Cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) if you haven't already. This is the engine to WDL.
2. Configure the AWS CLI for use with the rti-code AWS account. (configure with the secret key)
3. Clone the parent repository `biocloud_gwas_workflows` to your local machine.
   
```
cd /shared/
git clone https://github.com/RTIInternational/biocloud_gwas_workflows
```
   
4. Navigate to the `biocloud_gwas_workflows/meta_analysis/metal/ewas/local/` folder and edit the input file. (For instructions on this, see step three in the above AWS Batch section) 
* Note that all of your data/results must be local! So you will have to download your summary stats if they are on S3.
5. IMPORTANT STEP. Because you are running the workflow locally, you will have to edit the WDL scripts slightly. At the top of each WDL file, the location of each of the corresponding files is specified. This will have to be modified to reflect your working environment. Let me illustrate with an example.
   
<details>
<summary>full_ewas_meta.wdl details</summary>
                   
head of original file
```
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_preprocessing.wdl" as PREPROCESS
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_postprocessing.wdl" as POSTPROCESS
```            

head of edited file that will be ran locally
```
import "/shared/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS
import "/shared/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_preprocessing.wdl" as PREPROCESS
import "/shared/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_postprocessing.wdl" as POSTPROCESS
```      
Notice how I changed the paths to reflect where each file is in my local environment.
         
</details>
   
<details>
<summary>ewas_meta_postprocessing.wdl details</summary>
            
head of original file
```
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT
```
            
head of edited file that will be ran locally
```
import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS
import "/shared/jmarks/biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT
```
        
Notice how I changed the paths to reflect where each file is in my local environment.
</details>
   
<details>
<summary>ewas_meta_preprocessing.wdl details</summary>
         
head of original file
```
import "biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS       
```
            
head of edited file that will be ran locally
```
import "/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/ewas_meta_utils.wdl" as UTILS
```
        
Notice how I changed the paths to reflect where each file is in my local environment.
</details>
   
7. Execute the workflow. Make you have all of the data downloaded paths are correct to your input file.
   
```
cd /shared/biocloud_gwas_workflows/meta_analysis/metal/ewas/local/
java -jar ~/bin/cromwell/cromwell-54.jar run full_ewas_meta.wdl --inputs inputs.json
```
   
# Authors
   For any questions, comments, concerns, or bugs, send me an email or slack and I'll be happy to help.

   Jesse Marks (jmarks@rti.org)

