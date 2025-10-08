# Sex Check Workflow

Automated WDL workflow for performing a sex check.

## AWS Healthomics

### Start Run

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
** Parameter templates can be found [here](https://github.com/RTIInternational/biocloud_gwas_workflows/tree/master/genotype_qc_sub_wfs/sex_check/wdl_v1.1/config_templates/).
