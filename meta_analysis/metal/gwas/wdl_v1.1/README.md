# GWAS Meta-analysis Workflow: METAL

This repository provides a workflow for running a Genome-Wide Association Study (GWAS) meta-analysis using the [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) software tool.


## Workflow Overview

This workflow utilizes METAL to execute a GWAS meta-analysis. 

User needs GWAS summary statitics with the following columns (headers can be named whatever):
   - Variant ID
   - Chromosome
   - Position
   - Beta (effect_size)
   - SE (standard error of beta coefficient)
   - A1 (coded allele)
   - A2 (noncoded allele)
   - P-value
   - N

For each cohort, the results should be combined into one input file — i.e., merge the individual chromosome results.

## Running in AWS Healthomics

``` bash
# Set run parameters
host_dir="<Host directory to map to container>"
container_dir="<Container directory where host directory is mounted>"
aws_shared_credentials_file="<Path to AWS shared credential file>"
aws_profile="<AWS profile to use for credentials>"
charge_code="<Project charge code>"
workflow_id="<AWS Healthomics workflow ID>"
parameters="<JSON file containing run parameters>"
name="<Name of run>"
output_uri="<S3 location for run output>"
run_metadata_output_dir="<Dir where JSON file containing run metadata will be written>"
storage_capacity="<Size of storage for entire run>"

# Start run
docker run -ti \
    -v $host_dir:$container_dir \
    -e task=start_run \
    -e AWS_SHARED_CREDENTIALS_FILE="$aws_shared_credentials_file" \
    -e aws_profile="$aws_profile" \
    -e charge_code="$charge_code" \
    -e workflow_id=$workflow_id \
    -e parameters="$parameters" \
    -e name="$name" \
    -e output_uri="$output_uri" \
    -e run_metadata_output_dir="$run_metadata_output_dir" \
    -e storage_capacity=$storage_capacity \
    --rm rtibiocloud/healthomics_tools:v2.0_290ee15
```

## Config templates

- Same column names in all summary statistics files
    - [Template](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/meta_analysis/metal/gwas/wdl_v1.1/config_templates/metal_same_sum_stats_col_names_template.json)
    - [Example](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/meta_analysis/metal/gwas/wdl_v1.1/config_templates/metal_same_sum_stats_col_names_example.json)
- Different column names in summary statistics files
    - [Template](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/meta_analysis/metal/gwas/wdl_v1.1/config_templates/metal_different_sum_stats_col_names_template.json)
    - [Example](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/meta_analysis/metal/gwas/wdl_v1.1/config_templates/metal_different_sum_stats_col_names_example.json)

## Source

- [Main workflow](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/meta_analysis/metal/gwas/wdl_v1.1/metal_gwas_wf.wdl)
- WDL tools
    - [metal.wdl](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/metal/wdl_v1.1/metal.wdl)
    - [rti-tsv-utils.wdl](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/rti-tsv-utils/wdl_v1.1/rti-tsv-utils.wdl)
    - [tsv_utils.wdl](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/tsv_utils/wdl_v1.1/tsv_utils.wdl)
    - [generate_gwas_plots.wdl](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/generate_gwas_plots/wdl_v1.1/generate_gwas_plots.wdl)

## Authors
For any questions, please reach out to Jesse Marks (jmarks@rti.org) or Nathan Gaddis (ngaddis@rti.org).
