# Mahalanobis-Based Ancestry Prediction Workflow

## Introduction

## AWS Healthomics

### Create AWS Healthomics Workflow

``` bash
repo_dir="<REPO_DIR>"
docker run -ti \
    -v $repo_dir:$repo_dir \
    -e task=create_wf \
    -e aws_profile="<AWS_PROFILE>" \
    -e AWS_SHARED_CREDENTIALS_FILE="<AWS_SHARED_CREDENTIALS_FILE>" \
    -e repo_dir=$repo_dir \
    -e main=$repo_dir/ancestry/mahalanobis/wdl_v1.1/mahalanobis_ancestry_wf.wdl \
    -e name="mahalanobis_ancestry_wf" \
    -e description="Workflow for calculating the ancestral similarity of samples to reference populations based on Mahalanobis distance." \
    -e readme=$repo_dir/ancestry/mahalanobis/wdl_v1.1/README.md \
    -e engine=WDL \
    -e storage_capacity=2000 \
    --rm rtibiocloud/healthomics_tools:v2.0_7a0f3d4
```
Note: creation of workflow in AWS Healthomics only needs to be done before the first use and upon modifications to the workflow.

### Start Run
``` bash
docker run -ti \
    -v "<HOST_DIR>":"<CONTAINER_DIR>" \
    -e task=start_run \
    -e aws_profile="<AWS_PROFILE>" \
    -e AWS_SHARED_CREDENTIALS_FILE="<AWS_SHARED_CREDENTIALS_FILE>" \
    -e charge_code="<CHARGE_CODE>" \
    -e workflow_id="<WORKFLOW_ID>" \
    -e parameters="<PARAMETERS>" \
    -e name="<NAME>" \
    -e output_uri="<OUTPUT_URI>" \
    -e run_metadata_output_dir="<RUN_METADATA_OUTPUT_DIR>" \
    -e storage_capacity="<STORAGE_CAPACITY>" \
    --rm rtibiocloud/healthomics_tools:v2.0_7a0f3d4

```
