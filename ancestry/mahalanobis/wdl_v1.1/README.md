# Mahalanobis-Based Ancestry Prediction Workflow

## Introduction

## AWS Healthomics

### Start Run
``` bash
docker run -ti \
    -v "<HOST_DIR>":"<CONTAINER_DIR>" \
    -e task=start_run \
    -e aws_profile="<AWS_PROFILE>" \
    -e AWS_SHARED_CREDENTIALS_FILE="<AWS_SHARED_CREDENTIALS_FILE>" \
    -e charge_code="<CHARGE_CODE>" \
    -e workflow_id="<WORKFLOW_ID>" \
    -e workflow_version_name="<WORKFLOW_VERSION NAME>"
    -e name="<NAME>" \
    -e cache_id="<CACHE_ID>" \
    -e cache_behavior="<CACHE_BEHAVIOR>" \
    -e parameters="<PARAMETERS>" \
    -e output_uri="<OUTPUT_URI>" \
    -e run_metadata_output_dir="<RUN_METADATA_OUTPUT_DIR>" \
    -e workflow_type="<WORKFLOW_TYPE>" \
    -e priority="<PRIORITY>" \
    -e storage_type="<STORAGE_TYPE>" \
    -e storage_capacity="<STORAGE_CAPACITY>" \
    -e log_level="<LOG_LEVEL>" \
    -e retention_mode="<RETENTION_MODE>" \
    --rm rtibiocloud/healthomics_tools:v2.0_0076eb5
```

#### Parameters
| Parameter | Description | Type | Choices | Default Value | Required |
| --------- | ------ | ---- | ------- | ------------- | -------- |
| aws_profile | AWS profile to use for credentials | string  |  |  | Yes |
| AWS_SHARED_CREDENTIALS_FILE | Path to AWS shared credential file | string  |  |  | Yes |
| charge_code | RTI charge code | string |  |  | Yes |
| workflow_id | HealthOmics ID of workflow to run | string |  |  | Yes |
| workflow_version_name | Name of workflow version to run | String |  |  | No |
| name | Name to assign to run | string |  |  | Yes |
| cache_id | ID of cache to use for the run | string |  |  | No |
| cache_behavior | Cache behavior for the run | string | `CACHE_ON_FAILURE`, `CACHE_ALWAYS` |  | No |
| parameters | Path to JSON file containing run parameters | string |  |  | Yes |
| output_uri | S3 path for workflow output | string |  |  | Yes |
| run_metadata_output_dir | Directory to which run metadata will be output | string |  |  | Yes |
| workflow_type | Type of workflow to run | string |  `PRIVATE`, `READY2RUN` | `PRIVATE` | No |
| priority | Priority for run | integer | `1-100000` | `100` | No |
| storage_type | Storage type for run | string | `STATIC`, `DYNAMIC` | `STATIC` | No |
| storage_capacity | Storage capacity for run in GB if storage type = `STATIC` | integer | `1-10000` | `2000` | No |
| log_level | Log level for run | string | `OFF`, `FATAL`, `ERROR`, `ALL` | `ALL` | No |
| retention_mode | Retention mode for run | string | `RETAIN`, `REMOVE` | `RETAIN` | No |
