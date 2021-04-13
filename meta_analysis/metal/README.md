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

## Workflow overview
TBA

## Running the analysis
TBA




<br><br><br>
## To-do 
- [ ] Add 2df option
- [ ] Add option to apply [genomic control (GC)](https://en.wikipedia.org/wiki/Genomic_control).
- [ ] Add option to merge cohort specific GWAS results to the meta-analysis results.
- [ ] Add options to filter SNPs by imputation quality (R²) and minor allele frequency (MAF).
- [ ] Add option to opt out of removing singletons—the default is currently to remove singletons. 
- [ ] Create a tool to generate the input config file from an Excel file. 
