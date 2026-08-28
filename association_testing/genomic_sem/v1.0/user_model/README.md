# GenomicSEM User Model Workflow

## Overview

Fits a user-specified GenomicSEM structural equation model (SEM), defined via a lavaan model syntax file, against the genetic covariance structure of a set of traits. The workflow:

1. **Merges** per-chromosome summary statistics files within each study/trait (if multiple files are provided) into a single file prior to preprocessing using `tsv_append`.
2. **Preprocesses** each summary statistics file in parallel to standardize column names and extract rsIDs from variant identifiers.
3. **Munges** each preprocessed file in parallel using the reference SNP list.
4. Runs **LDSC** on the munged summary statistics to estimate the genetic covariance/correlation matrix.
5. Runs **usermodel()** to fit the user-specified lavaan model.

All memory-intensive tasks (munge, ldsc, usermodel) automatically size their runtime memory based on the size of their inputs unless an override is provided (see [Memory Calculation](#memory-calculation)).

## WDL File

[gsem_user_model.wdl](gsem_user_model.wdl)

## Workflow Steps

| Step | Task | Description |
| --- | --- | --- |
| 1 | `rti_tsv.tsv_append` (scattered, conditional) | Merges per-chromosome summary statistics files into a single file per study/trait when multiple files are provided. |
| 2 | `genomic_sem_preprocessing.genomic_sem_preprocessing` (scattered) | Standardizes column names and extracts rsIDs for each input sumstats file. |
| 3 | `gsem.gsem_munge` (scattered) | Munges each preprocessed sumstats file individually. |
| 4 | `gsem.gsem_ldsc` | Runs LDSC on all munged sumstats to produce the genetic covariance matrix. |
| 5 | `gsem.gsem_usermodel` | Fits the user-specified lavaan model. |

## Inputs

### Required

| Parameter | Type | Description |
| --- | --- | --- |
| `sumstats_files` | `Array[Array[File]]` | Input GWAS summary statistics files (either one genome sumstats file or per-chromosome sumstats files per trait/study). |
| `trait_names` | `Array[String]` | Trait labels corresponding to `sumstats_files`. |
| `sample_sizes` | `Array[Int]` | Per-trait sample sizes (N) for munging. |
| `ref_snp_list` | `File` | Reference SNP list (HM3) for munging. |
| `sumstats_columns` | `SUMSTATS_COLUMNS` | Per-file column-name mapping used for preprocessing (see [SUMSTATS_COLUMNS](#sumstats_columns-struct)). |
| `sample_prevs` | `Array[Float]` | Sample prevalences per trait for LDSC. |
| `population_prevs` | `Array[Float]` | Population prevalences per trait for LDSC. |
| `ld_dir` | `Directory` | Directory containing LD scores. |
| `wld_dir` | `Directory` | Directory containing LD score weights. |
| `model_lavaan` | `File` | Lavaan model syntax file for `usermodel()`. |

### Optional

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `estimation_method` | `String` | `"DWLS"` | Estimation method for `usermodel()`. |
| `info_filter` | `Float` | `0.9` | INFO/R2 threshold used by munge. |
| `maf_filter` | `Float` | `0.01` | MAF threshold used by munge. |
| `parallel` | `Boolean` | `false` | Enable parallel processing where supported. |
| `cores` | `Int` | `1` | Number of cores when parallel processing is enabled. |
| `no_overwrite` | `Boolean` | `false` | If true, do not overwrite existing munged outputs. |
| `cficalc` | `Boolean` | `false` | If true, compute CFI for usermodel output. |
| `std_lv` | `Boolean` | `false` | If true, use unit-variance identification for latent variables. |
| `imp_cov` | `Boolean` | `false` | If true, include implied/residual covariance output. |
| `fix_resid` | `Boolean` | `false` | If true, constrain residual variances above zero. |
| `toler` | `Float?` | — | Optional matrix inversion tolerance for usermodel. |
| `q_factor` | `Boolean` | `false` | If true, compute Q_Factor heterogeneity statistic. |
| `munge_out_dir` | `String` | `"munge_out"` | Output directory for munged summary statistics. |
| `preprocessing_out_dir` | `String` | `"preprocessing_out"` | Output directory for preprocessed summary statistics. |
| `ldsc_output_prefix` | `String` | `"ldsc/gsem_ldsc_output"` | Output prefix for LDSC outputs. |
| `usermodel_output_prefix` | `String` | `"usermodel/gsem_usermodel_output"` | Output prefix for usermodel output. |
| `ecr_repo` | `String?` | — | Optional ECR repository URI prefix used with each task's ECR image name. |
| `image_source` | `String` | `"docker"` | Container source selector: `docker` or `ecr`. |
| `tsv_append_cpu` | `Int` | `1` | CPU cores for the tsv_append task. |
| `tsv_append_mem_gb` | `Int` | `2` | Memory in GB for the tsv_append task. |
| `preprocessing_cpu` | `Int` | `1` | CPU cores for the preprocessing task. |
| `preprocessing_mem_gb` | `Int` | `4` | Memory in GB for the preprocessing task. |
| `munge_cpu` | `Int` | `1` | CPU cores for the munge task. |
| `munge_mem_gb` | `Int?` | — | Override for munge task memory (see [Memory Calculation](#memory-calculation)). |
| `ldsc_cpu` | `Int` | `1` | CPU cores for the ldsc task. |
| `ldsc_mem_gb` | `Int?` | — | Override for ldsc task memory (see [Memory Calculation](#memory-calculation)). |
| `usermodel_cpu` | `Int` | `1` | CPU cores for the usermodel task. |
| `usermodel_mem_gb` | `Int?` | — | Override for usermodel task memory (see [Memory Calculation](#memory-calculation)). |

### `SUMSTATS_COLUMNS` struct

Per-file column-name mapping used by the preprocessing step, indexed in the same order as `sumstats_files`:

```wdl
struct SUMSTATS_COLUMNS {
    Array[String] variant_id
    Array[String] effect_allele
    Array[String] non_effect_allele
    Array[String] effect
    Array[String] p
    Array[String?] z
    Array[String?] se
    Array[String?] n
    Array[String?] effect_allele_freq
    Array[String?] info
    Array[String?] direction
}
```

- `variant_id`, `effect_allele`, `non_effect_allele`, `effect`, and `p` are required for every file.
- `z`, `se`, `n`, `effect_allele_freq`, `info`, and `direction` are optional per-file (use `null` for files where a column is not available).

## Outputs

| Output | Type | Description |
| --- | --- | --- |
| `preprocessed_sumstats` | `Array[File]` | Preprocessed (column-standardized) summary statistics files. |
| `munged_sumstats` | `Array[File]` | Munged summary statistics files. |
| `ldsc_rds` | `File` | LDSC genetic covariance/correlation output (RDS). |
| `ldsc_log` | `File` | LDSC log file. |
| `usermodel_rds` | `File` | User model fit results (RDS). |

## Memory Calculation

For each of `munge`, `ldsc`, and `usermodel`, memory is automatically calculated as **1.5x the size of the task's inputs**, rounded up to `(2^n - 1)` GB, unless the corresponding `*_mem_gb` override is provided:

| Task | Memory basis when unset |
| --- | --- |
| munge | Per-file sumstats size + `ref_snp_list` size |
| ldsc | Munged sumstats size + `ld_dir` + `wld_dir` size (deduplicated if `ld_dir == wld_dir`) |
| usermodel | `ldsc_rds` + `model_lavaan` size |

## Notes

- `munge` and preprocessing run in parallel (scattered) across `sumstats_files`.
- `ecr_repo` and `image_source` are shared across every task call in the workflow; set `image_source = "ecr"` and provide `ecr_repo` to pull task images from ECR instead of Docker Hub.

## Example Inputs

See [gsem_user_model_inputs_min.json](gsem_user_model_inputs_min.json) for a minimal example and [gsem_user_model_inputs_full.json](gsem_user_model_inputs_full.json) for a full example with all optional parameters set.

## AWS HealthOmics Usage

The commands below use the [healthomics_tools](https://github.com/RTIInternational/biocloud_docker_tools/blob/master/healthomics_tools/v2.0/) [Docker image](https://hub.docker.com/repository/docker/rtibiocloud/healthomics_tools/tags) to register and run this workflow in AWS HealthOmics.

### Create Workflow

``` sh
# Define variables
BIOCLOUD_GWAS_WORKFLOWS_DIR="</path/to/biocloud_gwas_workflows>"
AWS_CREDENTIALS_DIR="</path/to/.aws>"
AWS_PROFILE="<aws_profile>"
VERSION="<vX.Y>"
HEALTHOMICS_TOOLS_DOCKER_IMAGE="rtibiocloud/healthomics_tools:<tag>"

# Create workflow
docker run -ti \
    -v "$BIOCLOUD_GWAS_WORKFLOWS_DIR":"$BIOCLOUD_GWAS_WORKFLOWS_DIR" \
    -v "$AWS_CREDENTIALS_DIR":"$AWS_CREDENTIALS_DIR" \
    -e task=create_wf \
    -e aws_profile="$AWS_PROFILE" \
    -e AWS_SHARED_CREDENTIALS_FILE="$AWS_CREDENTIALS_DIR/credentials" \
    -e repo_dir="$BIOCLOUD_GWAS_WORKFLOWS_DIR" \
    -e main="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/user_model/gsem_user_model.wdl" \
    -e name="gsem_user_model" \
    -e description="GenomicSEM user model workflow" \
    -e readme="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/user_model/README.md" \
    --rm $HEALTHOMICS_TOOLS_DOCKER_IMAGE
```

### Create Workflow Version

``` sh
# Define variables
BIOCLOUD_GWAS_WORKFLOWS_DIR="</path/to/biocloud_gwas_workflows>"
AWS_CREDENTIALS_DIR="</path/to/.aws>"
AWS_PROFILE="<aws_profile>"
WORKFLOW_ID="<workflow_id>"
VERSION="<vX.Y>"
HEALTHOMICS_TOOLS_DOCKER_IMAGE="rtibiocloud/healthomics_tools:<tag>"

# Create workflow version
docker run -ti \
    -v "$BIOCLOUD_GWAS_WORKFLOWS_DIR":"$BIOCLOUD_GWAS_WORKFLOWS_DIR" \
    -v "$AWS_CREDENTIALS_DIR":"$AWS_CREDENTIALS_DIR" \
    -e task=create_wf_version \
    -e aws_profile="$AWS_PROFILE" \
    -e AWS_SHARED_CREDENTIALS_FILE="$AWS_CREDENTIALS_DIR/credentials" \
    -e workflow_id="$WORKFLOW_ID" \
    -e repo_dir="$BIOCLOUD_GWAS_WORKFLOWS_DIR" \
    -e main="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/user_model/gsem_user_model.wdl" \
    -e name="gsem_user_model_$VERSION" \
    -e description="GenomicSEM user model workflow" \
    -e readme="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/user_model/README.md" \
    --rm $HEALTHOMICS_TOOLS_DOCKER_IMAGE
```

### Start Run

``` sh
# Define variables
DATA_DIR="</path/to/data/dir>"
AWS_CREDENTIALS_DIR="</host/path/to/.aws>"
CHARGE_CODE="<charge_code>"
AWS_PROFILE="<aws_profile>"
WORKFLOW_ID="<workflow_id>"
WORKFLOW_VERSION_NAME="<workflow_version_name>"
CACHE_ID="<cache_id>"
CACHE_BEHAVIOR="<CACHE_ON_FAILURE|CACHE_ALWAYS>"
RUN_NAME="<run_name>"
JSON_INPUTS_PATH="$DATA_DIR/<path/to/json/inputs>"
OUTPUT_URI="</s3/path/for/workflow/output>"
RUN_METADATA_OUTPUT_DIR="$DATA_DIR/<path/to/run_metadata_output_dir>"
HEALTHOMICS_TOOLS_DOCKER_IMAGE="rtibiocloud/healthomics_tools:<tag>"

# Start run
docker run -ti \
    -u $(id -u):$(id -g) \
    -v "$DATA_DIR":"$DATA_DIR" \
    -v "$AWS_CREDENTIALS_DIR":"$AWS_CREDENTIALS_DIR" \
    -e task=start_run \
    -e charge_code="$CHARGE_CODE" \
    -e aws_profile="$AWS_PROFILE" \
    -e AWS_SHARED_CREDENTIALS_FILE="$AWS_CREDENTIALS_DIR/credentials" \
    -e workflow_id="$WORKFLOW_ID" \
    -e workflow_version_name="$WORKFLOW_VERSION_NAME" \
    -e cache_id="$CACHE_ID" \
    -e cache_behavior="$CACHE_BEHAVIOR" \
    -e name="$RUN_NAME" \
    -e parameters="$JSON_INPUTS_PATH" \
    -e output_uri="$OUTPUT_URI" \
    -e run_metadata_output_dir="$RUN_METADATA_OUTPUT_DIR" \
    --rm $HEALTHOMICS_TOOLS_DOCKER_IMAGE
```

#### Notes

- See the [healthomics_tools README](https://github.com/RTIInternational/biocloud_docker_tools/blob/master/healthomics_tools/v2.0/README.md) for the full list of optional parameters for each command (for example, `storage_capacity`, `cache_id`, `priority`, `role_arn`).
- To use the `cache_id`/`cache_behavior` parameters, a run cache must already exist; create one first with [Create Run Cache](https://github.com/RTIInternational/biocloud_docker_tools/blob/master/healthomics_tools/v2.0/README.md#create-run-cache) and pass its `id` as `cache_id`.
- See [Docker Hub](https://hub.docker.com/repository/docker/rtibiocloud/healthomics_tools/tags) for the most recent HealthOmics Tools image tag
