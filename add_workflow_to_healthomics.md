# HealthOmics Workflow Setup

## Docker images
### Create ECR images for all Docker images used in workflow if they don't already exist
Example:
```bash
# aws_account_id=404545384114  # Code
aws_account_id=515876044319  #IOmics
aws_region=us-east-1
aws_profile=<AWS_PROFILE>
repository_tag=<REPO_TAG>
project_number=<PROJECT_NUMBER>

# Pull image to be transferred to local system
docker pull $repository_tag

# Tag image
docker tag $repository_tag $aws_account_id.dkr.ecr.us-east-1.amazonaws.com/$repository_tag

# Create repository in ECR for image if it doesn't exist
repository_name=$(echo $repository_tag | perl -pe 's/:.+//;')
aws ecr create-repository \
  --repository-name $repository_name \
  --tags '[{"Key": "project-number","Value": "'$project_number'"}]' \
  --profile $aws_profile

# Set permissions on ECR repo
policy_json=$(cat <<'EOF'
{
  "Version" : "2012-10-17",
  "Statement" : [ {
    "Sid" : "omics workflow",
    "Effect" : "Allow",
    "Principal" : {
      "Service" : "omics.amazonaws.com"
    },
    "Action" : [ "ecr:BatchCheckLayerAvailability", "ecr:BatchGetImage", "ecr:GetDownloadUrlForLayer" ]
  } ]
}
EOF
)
aws ecr set-repository-policy \
    --repository-name $repository_name \
    --policy-text "$policy_json" \
    --profile $aws_profile


# Authenticate your Docker client to the Amazon ECR registry
aws ecr get-login-password --region $aws_region --profile $aws_profile | docker login --username AWS --password-stdin $aws_account_id.dkr.ecr.$aws_region.amazonaws.com

# Push image
docker push $aws_account_id.dkr.ecr.us-east-1.amazonaws.com/$repository_tag
```

## WDL
### Update version
Add `version 1.1` at beginning of WDL file

### Update container specifications
Example:
```wdl
String docker_image = "rtibiocloud/assign_ancestry_mahalanobis:v1_5a4fd59"
String ecr_image = "rtibiocloud/assign_ancestry_mahalanobis:v1_5a4fd59"
String? ecr_repo
String image_source = "docker"
String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
```

### Change WDL variables from `${}` syntax to `~{}` syntax
Old syntax:
``` wdl
${chr}
```
New syntax:
``` wdl
~{chr}
```

### Change syntax of `sep` function
Old syntax:
``` wdl
~{sep="," input_files}
```
New syntax:
``` wdl
~{sep(',', input_files)}
```

### Change ternary operator
Old syntax:
``` wdl
${true='--snps-only' false='' snps_only}
```
New syntax:
``` wdl
 ~{if snps_only then '--snps-only' else ''}
```

### Remove paths from imports
Old syntax:
```wdl
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT
```
New syntax:
``` wdl
import "convert_variant_ids.wdl" as IDCONVERT
```

### Remove `maxRetries` runtime parameter
Old syntax:
```wdl
Int max_retries = 3

runtime{
    docker: container_image
    cpu: cpu
    memory: "${mem_gb} GB"
    maxRetries: max_retries
}
```
New syntax:
```wdl
runtime{
    docker: container_image
    cpu: cpu
    memory: "${mem_gb} GB"
}
```

### Enclose input parameters in input block
Old syntax:
```wdl
File summary_stats
String col_id
String col_chromosome
String col_position
String col_p
String output_basename
```
New syntax:
```wdl
input {
    File summary_stats
    String col_id
    String col_chromosome
    String col_position
    String col_p
    String output_basename
}
```

## Create a dependencies file for the workflow and sub-workflows
The dependencies file is a json file in the same directory as the workflow/sub-workflow. The name of the file should be `<WORKFLOW_NAME>_dependencies.json`.
Example dependencies file:
```json
{
    "workflows": [
        "genotype_qc_sub_wfs/ld_pruning/wdl_v1.1/ld_prune_wf.wdl"

    ],
    "wdl_tools": [
        "biocloud_wdl_tools/assign_ancestry_mahalanobis/wdl_v1.1/assign_ancestry_mahalanobis.wdl",
        "biocloud_wdl_tools/convert_variant_ids/wdl_v1.1/convert_variant_ids.wdl",
        "biocloud_wdl_tools/plink/wdl_v1.1/plink.wdl",
        "biocloud_wdl_tools/smartpca/wdl_v1.1/smartpca.wdl",
        "biocloud_wdl_tools/utils/wdl_v1.1/utils.wdl"
    ],
    "structs": []
}
```

## Create a parameters file for the workflow 
The parameters file is a json file in the same directory as the workflow/sub-workflow. The name of the file should be `<WORKFLOW_NAME>_parameters.json`.
Example parameters file:
```json
{
	"dataset_bed": {
		"description": "Dataset bed file"
	},
	"dataset_bim": {
		"description": "Dataset bim file"
	},
	"dataset_fam": {
		"description": "Dataset fam file"
	},
    	"ancestry_pop_type": {
		"description": "Population type (superpop or pop)",
        "optional": true
	},
	"ancestries_to_include": {
		"description": "Array of ancestries to include",
        "optional": true
	}
}
```

## Create a README file for the workflow
Create a `README.md` file for the workflow that describes the workflow and provides an example command for running the workflow.

## Add workflow to HealthOmics
Example:
```bash
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
    -e main=$repo_dir/ancestry/mahalanobis/wdl_v1.1/mahalanobis_ancestry_wf.wdl \
    -e parameter_template=$repo_dir/ancestry/mahalanobis/wdl_v1.1/mahalanobis_ancestry_wf_parameters.json \
    -e name="mahalanobis_ancestry_wf" \
    -e description="Workflow for calculating the ancestral similarity of samples to reference populations based on Mahalanobis distance." \
    -e engine=WDL \
    -e storage_capacity=2000 \
    --rm rtibiocloud/healthomics_tools:v1.0_f456174
```
