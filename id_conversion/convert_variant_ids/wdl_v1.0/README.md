### Create Workflow in AWS HealthOmics

``` bash
aws_profile=iomics
aws_access_key_id=$(perl -ane 'BEGIN {$take = 0;} if ($F[0] =~ /'$aws_profile'/) { $take = 1; } if ($take && $F[0] =~ /aws_access_key_id/) { print $F[2]; $take = 0; }' ~/.aws/credentials)
aws_secret_access_key=$(perl -ane 'BEGIN {$take = 0;} if ($F[0] =~ /'$aws_profile'/) { $take = 1; } if ($take && $F[0] =~ /aws_secret_access_key/) { print $F[2]; $take = 0; }' ~/.aws/credentials)
aws_region_name="us-east-1"
repo_dir="/home/ngaddis/git/biocloud_gwas_workflows"
main="id_conversion/convert_variant_ids/wdl_v1.0/convert_variant_ids_wf.wdl"
name="convert_variant_ids"
description="Workflow for converting variant IDs to standard format"
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