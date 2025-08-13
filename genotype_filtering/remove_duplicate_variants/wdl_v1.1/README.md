# Remove Duplicate Variants Workflow

## Introduction

## AWS Healthomics

### Create AWS Healthomics Workflow

``` bash
aws_profile=code
aws_access_key_id=$(perl -ane 'BEGIN {$take = 0;} if ($F[0] =~ /'$aws_profile'/) { $take = 1; } if ($take && $F[0] =~ /aws_access_key_id/) { print $F[2]; $take = 0; }' ~/.aws/credentials)
aws_secret_access_key=$(perl -ane 'BEGIN {$take = 0;} if ($F[0] =~ /'$aws_profile'/) { $take = 1; } if ($take && $F[0] =~ /aws_secret_access_key/) { print $F[2]; $take = 0; }' ~/.aws/credentials)
aws_region_name="us-east-1"
repo_dir="/rti-01/ngaddis/git/biocloud_gwas_workflows"
main="genotype_filtering/remove_duplicate_variants/wdl_v1.1/remove_duplicate_variants_wf.wdl"
name="remove_duplicate_variants"
description="Workflow for removing duplicate variants from genotypes in plink format"
docker run -ti \
    -v $repo_dir:$repo_dir \
    -e task=create_wf \
    -e aws_access_key_id="$aws_access_key_id" \
    -e aws_secret_access_key="$aws_secret_access_key" \
    -e aws_region_name="$aws_region_name"\
    -e repo_dir="$repo_dir" \
    -e main="$main" \
    -e name="$name" \
    -e description="$description" \
    -e engine="WDL" \
    -e storage_capacity=2000 \
    --rm rtibiocloud/healthomics_tools:v1.0_f456174
```
Note: creation of workflow in AWS Healthomics only needs to be done before the first use and upon modifications to the workflow.


### Start Run
``` bash
host_dir="<Host directory to map to container>"
container_dir="<Container directory where host directory is mounted>"
charge_code="<Project charge code>"
aws_access_key_id="<AWS access key ID>"
aws_secret_access_key="<AWS secret access key>"
aws_region_name="<AWS region, e.g., us-east-1>"
workflow_id="<ID of workflow to run>"
parameters="<JSON file containing run parameters>"
name="<Name of run>"
output_uri="<S3 location for run output>"
run_metadata_output_dir="<Dir where JSON file containing run metadata will be written>"
storage_capacity="<Size of storage for entire run>"
docker run -ti \
    -v $host_dir:$container_dir \
    -e task=start_run \
    -e charge_code=$charge_code \
    -e aws_access_key_id=$aws_access_key_id \
    -e aws_secret_access_key=$aws_secret_access_key \
    -e aws_region_name=$aws_region_name \
    -e workflow_id=$workflow_id \
    -e parameters=$parameters \
    -e name=$name \
    -e output_uri=$output_uri \
    -e run_metadata_output_dir=$run_metadata_output_dir \
    -e storage_capacity=$storage_capacity \
    --rm rtibiocloud/healthomics_tools:v1.0_f456174

```

