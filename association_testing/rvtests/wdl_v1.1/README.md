# Rvtests GWAS Workflow

Automated WDL workflow for the GWAS software [Rvtests (Rare variant tests)](http://zhanxw.github.io/rvtests/).
[Source code](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/association_testing/rvtests/wdl_v1.1/)


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

- Autosomes
    - Unrelated
        - Continuous
            - [CODE](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_unrelated_continuous_code.json)
            - [IOMICS](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_unrelated_continuous_iomics.json)
        - Logistic
            - [CODE](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_unrelated_logistic_code.json)
            - [IOMICS](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_unrelated_logistic_iomics.json)
    - Related
        - Continuous
            - [CODE](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_related_continuous_code.json)
            - [IOMICS](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_related_continuous_iomics.json)
        - Logistic
            - [CODE](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_related_logistic_code.json)
            - [IOMICS](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/config_templates/rvtests_autosomes_related_logistic_iomics.json)

## Source

- [Main workflow](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/association_testing/rvtests/wdl_v1.1/rvtests_gwas_wf.wdl)
- Sub-workflows
    - [rvtests_gwas_chr_wf](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/association_testing/rvtests/wdl_v1.1/rvtests_gwas_chr_wf.wdl)
    - [summarize_gwas_wf](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/helper_workflows/summarize_gwas/wdl_v1.1/summarize_gwas_wf.wdl)
    - [generate_kinship_matrix_wf](https://github.com/RTIInternational/biocloud_gwas_workflows/blob/master/helper_workflows/generate_kinship_matrix/wdl_v1.1/generate_kinship_matrix_wf.wdl)
- WDL tools
    - [rti-tsv-utils](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/rti-tsv-utils/wdl_v1.1/rti-tsv-utils.wdl)
    - [convert_variant_ids.wdl](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/convert_variant_ids/wdl_v1.1/convert_variant_ids.wdl)
    - [utils.wdl](https://github.com/RTIInternational/biocloud_wdl_tools/blob/master/utils/wdl_v1.1/utils.wdl)
