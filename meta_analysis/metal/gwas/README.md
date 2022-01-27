# GWAS meta-analysis workflow: METAL
Documentation for running a GWAS meta-analysis using the [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) software tool.

## Pre-requisites
* Unix-based operating system (Linux or OSx. Sorry Windows folks.)
* Java v1.8 and higher (download here)
* Docker
* AWS CLI

1. Install [Cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) if you haven't already
2. [Configure the AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-configure.html) for use with CODE AWS group.
    * Configure with the secret key associated with the CODE AWS account
4. Clone local copy of biocloud_gwas_workflows to download this repository to your local machine
```
git clone https://github.com/RTIInternational/biocloud_gwas_workflows
```
Now you're ready to get started running the workflow.

More detailed running instructions for WDL/Cromwell workflows are found in the tutorial for the metaxcan pipeline [here](https://github.com/RTIInternational/metaxcan-pipeline).




<br><br>

## Workflow overview
User needs GWAS summary statitics with the following columns (headers can be named whatever):
* rsID
* chromosome
* position
* Beta (effect_size)
* se_beta (standard error of beta coefficient)
* A1 (coded allele)
* A2 (noncoded allele)
* P-value

For each cohort the results should be combined into one input file—i.e. merge the individual chromosome results.

## Running the analysis
See the Microsoft Teams wiki [WDL/Cromwell Cheat Sheet](https://teams.microsoft.com/l/entity/com.microsoft.teamspace.tab.wiki/tab::61aecad5-13fa-4bde-adce-ba3b16950439?context=%7B%22subEntityId%22%3A%22%7B%5C%22pageId%5C%22%3A18%2C%5C%22origin%5C%22%3A2%7D%22%2C%22channelId%22%3A%2219%3Af42632e48b7c4b9e9f362afa1e4e1957%40thread.tacv2%22%7D&tenantId=2ffc2ede-4d44-4994-8082-487341fa43fb) for a reference. 

### Input file

| chr | pos | name | A1 | A2 |
|-----|-----|------|----|----|


All column number entries are 1-based. For example, the above portion of a summary stats header would have
* metal_gwas_meta_analysis_wf:variant_id_column: 3
* metal_gwas_meta_analysis_wf:chromosome_column: 1
* metal_gwas_meta_analysis_wf:pos_column: 2

remove_singletons can be either "true" or "false" (don't forget the quotations!). Setting this to true removes any variants that were present in only one of the input files. 

chromosomes_to_keep is a list of integers denoting the which chromosomes to run the analysis on. This allows you to run a meta based on a subset of chromosomes if desired.


<br><br><br>
## To-do 
- [ ] Add 2df option
- [ ] Add option to apply [genomic control (GC)](https://en.wikipedia.org/wiki/Genomic_control).


# Authors
For any questions, comments, concerns, or bugs, send me an email or slack and I'll be happy to help.

Jesse Marks (jmarks@rti.org)
