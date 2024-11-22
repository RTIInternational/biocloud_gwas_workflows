# Pulmonary Longitudinal GxE Workflow

## Introduction

## AWS Healthomics

### Create AWS Healthomics Workflow

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
    -e main=$repo_dir/association_testing/pulmonary_longitudinal_gxe/wdl_v1.1/pulmonary_longitudinal_gxe_wf.wdl \
    -e parameter_template=$repo_dir/association_testing/pulmonary_longitudinal_gxe/wdl_v1.1/pulmonary_longitudinal_gxe_wf_parameters.json \
    -e name="pulmonary_longitudinal_gxe" \
    -e description="Workflow for performing longitudinal GxE GWAS of lung function and omega3." \
    -e engine=WDL \
    -e storage_capacity=2000 \
    --rm rtibiocloud/healthomics_tools:v1.0_f456174
```
Note: creation of workflow in AWS Healthomics only needs to be done before the first use and upon modifications to the workflow.

### Start Run
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
