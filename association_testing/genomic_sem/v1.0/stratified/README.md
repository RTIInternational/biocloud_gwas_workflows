# GenomicSEM Stratified LDSC Enrichment Workflow

## Overview

Runs stratified (partitioned) LD score regression across functional annotation categories and estimates enrichment of a user-specified GenomicSEM model's parameters within those categories. The workflow:

1. **Merges** per-chromosome summary statistics files within each study/trait (if multiple files are provided) into a single file prior to preprocessing using `tsv_append`.
2. **Preprocesses** each summary statistics file in parallel to standardize column names and extract rsIDs from variant identifiers.
3. **Munges** each preprocessed file in parallel using the reference SNP list.
4. Runs **stratified LDSC (s_ldsc)** on the munged summary statistics using partitioned LD scores.
5. Runs **enrich()** to estimate enrichment of the specified lavaan model parameters across annotation categories.

All memory-intensive tasks (munge, s_ldsc, enrich) automatically size their runtime memory based on the size of their inputs unless an override is provided (see [Memory Calculation](#memory-calculation)).

## WDL File

[gsem_stratified.wdl](gsem_stratified.wdl)

## Workflow Steps

| Step | Task | Description |
| --- | --- | --- |
| 1 | `rti_tsv.tsv_append` (scattered, conditional) | Merges per-chromosome summary statistics files into a single file per study/trait when multiple files are provided. |
| 2 | `genomic_sem_preprocessing.genomic_sem_preprocessing` (scattered) | Standardizes column names and extracts rsIDs for each input sumstats file. |
| 3 | `gsem.gsem_munge` (scattered) | Munges each preprocessed sumstats file individually. |
| 4 | `gsem.gsem_s_ldsc` | Runs stratified LDSC using partitioned LD scores, weights, and allele frequencies. |
| 5 | `gsem.gsem_enrich` | Estimates enrichment of the specified lavaan model parameters. |

## Inputs

### Required

| Parameter | Type | Description |
| --- | --- | --- |
| `sumstats_files` | `Array[Array[File]]` | Input GWAS summary statistics files (either one genome sumstats file or per-chromosome sumstats files per trait/study). |
| `trait_names` | `Array[String]` | Trait labels corresponding to `sumstats_files`. |
| `sample_sizes` | `Array[Int]` | Per-trait sample sizes (N) for munging. |
| `ref_snp_list` | `File` | Reference SNP list (HM3) for munging. |
| `sumstats_columns` | `SUMSTATS_COLUMNS` | Per-file column-name mapping used for preprocessing (see [SUMSTATS_COLUMNS](#sumstats_columns-struct)). |
| `sample_prevs` | `Array[Float]` | Sample prevalences per trait for stratified LDSC. |
| `population_prevs` | `Array[Float]` | Population prevalences per trait for stratified LDSC. |
| `ld_dir` | `Directory` | Directory containing partitioned LD scores. |
| `wld_dir` | `Directory` | Directory containing LD score weights. |
| `frq_dir` | `Directory` | Directory containing allele frequency files for s_ldsc. |
| `model_lavaan` | `File` | Lavaan model syntax file for `enrich()`. |
| `params` | `File` | File listing target lavaan parameters for enrichment estimation. |

### Optional

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `fixparam` | `File?` | — | Optional file listing parameters to fix during enrichment estimation. |
| `info_filter` | `Float` | `0.9` | INFO/R2 threshold used by munge. |
| `maf_filter` | `Float` | `0.01` | MAF threshold used by munge. |
| `parallel` | `Boolean` | `false` | Enable parallel processing where supported. |
| `cores` | `Int` | `1` | Number of cores when parallel processing is enabled. |
| `no_overwrite` | `Boolean` | `false` | If true, do not overwrite existing munged outputs. |
| `n_blocks` | `Int` | `200` | Number of jackknife blocks for s_ldsc standard error estimation. |
| `include_cont` | `Boolean` | `false` | If true, include continuous annotations in s_ldsc. |
| `fix` | `String` | `"regressions"` | Parameter class to fix in the enrichment model: `regressions`, `variances`, or `covariances`. |
| `std_lv` | `Boolean` | `false` | If true, use unit-variance identification for latent variables. |
| `not_rm_flank` | `Boolean` | `false` | If true, keep flanking-window and continuous annotations in enrich output. |
| `tau` | `Boolean` | `false` | If true, use tau matrices instead of zero-order matrices in enrich. |
| `not_base` | `Boolean` | `false` | If true, exclude baseline model estimates from enrich output. |
| `toler` | `Float?` | — | Optional matrix inversion tolerance for enrich. |
| `munge_out_dir` | `String` | `"munge_out"` | Output directory for munged summary statistics. |
| `preprocessing_out_dir` | `String` | `"preprocessing_out"` | Output directory for preprocessed summary statistics. |
| `s_ldsc_output_prefix` | `String` | `"s_ldsc/gsem_s_ldsc_output"` | Output prefix for stratified LDSC outputs. |
| `enrich_output_prefix` | `String` | `"enrich/gsem_enrich_output"` | Output prefix for enrich output. |
| `ecr_repo` | `String?` | — | Optional ECR repository URI prefix used with each task's ECR image name. |
| `image_source` | `String` | `"docker"` | Container source selector: `docker` or `ecr`. |
| `tsv_append_cpu` | `Int` | `1` | CPU cores for the tsv_append task. |
| `tsv_append_mem_gb` | `Int` | `2` | Memory in GB for the tsv_append task. |
| `preprocessing_cpu` | `Int` | `1` | CPU cores for the preprocessing task. |
| `preprocessing_mem_gb` | `Int` | `4` | Memory in GB for the preprocessing task. |
| `munge_cpu` | `Int` | `1` | CPU cores for the munge task. |
| `munge_mem_gb` | `Int?` | — | Override for munge task memory (see [Memory Calculation](#memory-calculation)). |
| `s_ldsc_cpu` | `Int` | `1` | CPU cores for the s_ldsc task. |
| `s_ldsc_mem_gb` | `Int?` | — | Override for s_ldsc task memory (see [Memory Calculation](#memory-calculation)). |
| `enrich_cpu` | `Int` | `1` | CPU cores for the enrich task. |
| `enrich_mem_gb` | `Int?` | — | Override for enrich task memory (see [Memory Calculation](#memory-calculation)). |

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
| `s_ldsc_rds` | `File` | Stratified LDSC output (RDS). |
| `s_ldsc_log` | `File` | Stratified LDSC log file. |
| `enrich_rds` | `File` | Enrichment estimation results (RDS). |

## Memory Calculation

For each of `munge`, `s_ldsc`, and `enrich`, memory is automatically calculated as **1.5x the size of the task's inputs**, rounded up to `(2^n - 1)` GB, unless the corresponding `*_mem_gb` override is provided:

| Task | Memory basis when unset |
| --- | --- |
| munge | Per-file sumstats size + `ref_snp_list` size |
| s_ldsc | Munged sumstats size + `ld_dir` + `wld_dir` + `frq_dir` size (deduplicated across any directories that are identical) |
| enrich | `s_ldsc_rds` + `model_lavaan` + `params` size |

## Notes

- `munge` and preprocessing run in parallel (scattered) across `sumstats_files`.
- `ecr_repo` and `image_source` are shared across every task call in the workflow; set `image_source = "ecr"` and provide `ecr_repo` to pull task images from ECR instead of Docker Hub.

## Example Inputs

See [gsem_stratified_inputs_min.json](gsem_stratified_inputs_min.json) for a minimal example and [gsem_stratified_inputs_full.json](gsem_stratified_inputs_full.json) for a full example with all optional parameters set.

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
    -e main="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/stratified/gsem_stratified.wdl" \
    -e name="gsem_stratified" \
    -e description="GenomicSEM stratified LDSC enrichment workflow" \
    -e readme="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/stratified/README.md" \
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
    -e main="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/stratified/gsem_stratified.wdl" \
    -e name="gsem_stratified_$VERSION" \
    -e description="GenomicSEM stratified LDSC enrichment workflow" \
    -e readme="$BIOCLOUD_GWAS_WORKFLOWS_DIR/association_testing/genomic_sem/$VERSION/stratified/README.md" \
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
