version development

import "genomic_sem.wdl" as gsem
import "utils.wdl" as utils
import "genomic_sem_preprocessing.wdl" as preprocessing
import "rti-tsv-utils.wdl" as rti_tsv

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

workflow gsem_common_factor_model {
  input {
    Array[Array[File]] sumstats_files
    Array[String] trait_names
    Array[Int] sample_sizes
    File ref_snp_list

    SUMSTATS_COLUMNS sumstats_columns

    Array[Float] sample_prevs
    Array[Float] population_prevs
    Directory ld_dir
    Directory wld_dir

    String estimation_method = "DWLS"

    Float info_filter = 0.9
    Float maf_filter = 0.01
    Boolean parallel = false
    Int cores = 1
    Boolean no_overwrite = false

    String munge_out_dir = "munge_out"
    String ldsc_output_prefix = "ldsc/gsem_ldsc_output"
    String commonfactor_output_prefix = "commonfactor/gsem_commonfactor_output"

    String? ecr_repo
    String image_source = "docker"
    Int tsv_append_cpu = 1
    Int tsv_append_mem_gb = 2
    String preprocessing_out_dir = "preprocessing_out"
    Int preprocessing_cpu = 1
    Int preprocessing_mem_gb = 4
    Int munge_cpu = 1
    Int? munge_mem_gb
    Int ldsc_cpu = 1
    Int? ldsc_mem_gb
    Int commonfactor_cpu = 1
    Int? commonfactor_mem_gb
}

  parameter_meta {
    sumstats_files: "Input GWAS summary statistics files (either one genome sumstats file or per-chromosome sumstats files per trait/study)."
    trait_names: "Trait labels corresponding to sumstats_files."
    sample_sizes: "Per-trait sample sizes (N) for munging."
    ref_snp_list: "Reference SNP list (HM3) for munging."
    sample_prevs: "Sample prevalences per trait for LDSC."
    population_prevs: "Population prevalences per trait for LDSC."
    ld_dir: "Directory containing LD scores."
    wld_dir: "Directory containing LD score weights."
    estimation_method: "Estimation method for commonfactor() (for example, DWLS)."
    info_filter: "INFO/R2 threshold used by gsem_munge."
    maf_filter: "MAF threshold used by gsem_munge."
    parallel: "Enable parallel processing where supported."
    cores: "Number of cores when parallel processing is enabled."
    no_overwrite: "If true, do not overwrite existing munged outputs."
    munge_out_dir: "Output directory for munged summary statistics."
    ldsc_output_prefix: "Output prefix for LDSC outputs."
    commonfactor_output_prefix: "Output prefix for common factor model output."
    ecr_repo: "Optional ECR repository URI prefix used with the task ecr image name."
    image_source: "Container source selector: docker or ecr."
    tsv_append_cpu: "Requested CPU cores for tsv_append task runtime."
    tsv_append_mem_gb: "Memory in GB for tsv_append task runtime."
    preprocessing_out_dir: "Output directory for preprocessed summary statistics."
    sumstats_columns: "Column mapping for preprocessing."
    preprocessing_cpu: "Requested CPU cores for preprocessing task runtime."
    preprocessing_mem_gb: "Memory in GB for preprocessing task runtime."
    munge_cpu: "Requested CPU cores for munge task runtime."
    munge_mem_gb: "Optional override for munge task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size), rounded to (2^n - 1) GB."
    ldsc_cpu: "Requested CPU cores for ldsc task runtime."
    ldsc_mem_gb: "Optional override for ldsc task runtime memory; if unset, memory is computed as 1.5x(munged sumstats + ld_dir + wld_dir sizes), rounded to (2^n - 1) GB."
    commonfactor_cpu: "Requested CPU cores for commonfactor task runtime."
    commonfactor_mem_gb: "Optional override for commonfactor task runtime memory; if unset, memory is computed as 1.5x(ldsc_rds size), rounded to (2^n - 1) GB."
}

  Int num_traits = length(sumstats_files)
  Boolean input_lengths_match = (
    length(trait_names) == num_traits &&
    length(sample_sizes) == num_traits &&
    length(sample_prevs) == num_traits &&
    length(population_prevs) == num_traits &&
    length(sumstats_columns.variant_id) == num_traits &&
    length(sumstats_columns.effect_allele) == num_traits &&
    length(sumstats_columns.non_effect_allele) == num_traits &&
    length(sumstats_columns.effect) == num_traits &&
    length(sumstats_columns.p) == num_traits &&
    length(sumstats_columns.z) == num_traits &&
    length(sumstats_columns.se) == num_traits &&
    length(sumstats_columns.n) == num_traits &&
    length(sumstats_columns.effect_allele_freq) == num_traits &&
    length(sumstats_columns.info) == num_traits &&
    length(sumstats_columns.direction) == num_traits
  )

  if (!input_lengths_match) {
    call utils.raise_error as input_length_mismatch_error {
      input:
        msg = "Workflow input error: sumstats_files, trait_names, sample_sizes, sample_prevs, population_prevs, and all arrays within sumstats_columns must have the same length.",
        ecr_repo = ecr_repo,
        image_source = image_source
    }
  }

  scatter (idx in range(num_traits)) {
    if (length(sumstats_files[idx]) > 1) {
      call rti_tsv.tsv_append as merge_sumstats {
        input:
          input_files = sumstats_files[idx],
          output_prefix = "~{preprocessing_out_dir}/merged_sumstats_~{idx}",
          ecr_repo = ecr_repo,
          image_source = image_source,
          cpu = tsv_append_cpu,
          mem_gb = tsv_append_mem_gb
      }
    }

    File trait_sumstats = select_first([merge_sumstats.out_tsv, sumstats_files[idx][0]])

    call preprocessing.genomic_sem_preprocessing as preprocess {
      input:
        sumstats_file = trait_sumstats,
        col_variant_id = sumstats_columns.variant_id[idx],
        col_effect_allele = sumstats_columns.effect_allele[idx],
        col_non_effect_allele = sumstats_columns.non_effect_allele[idx],
        col_effect = sumstats_columns.effect[idx],
        col_p = sumstats_columns.p[idx],
        col_z = sumstats_columns.z[idx],
        col_se = sumstats_columns.se[idx],
        col_n = sumstats_columns.n[idx],
        col_effect_allele_freq = sumstats_columns.effect_allele_freq[idx],
        col_info = sumstats_columns.info[idx],
        col_direction = sumstats_columns.direction[idx],
        out_file = "~{preprocessing_out_dir}/preprocessed_~{idx}",
        docker_image = "genomic-sem-preprocessing:v1",
        ecr_image = "genomic-sem-preprocessing:v1",
        ecr_repo = ecr_repo,
        image_source = image_source,
        cpu = preprocessing_cpu,
        mem_gb = preprocessing_mem_gb
    }
  }

  scatter (idx in range(num_traits)) {
    if (!defined(munge_mem_gb)) {
      call utils.get_total_file_size as munge_input_size {
        input:
          input_files = [preprocess.processed_sumstats[idx], ref_snp_list]
      }

      call utils.round_power_of_two_minus_one as calc_munge_memory {
        input:
          input_value = munge_input_size.total_file_size_gb * 1.5
      }
      Int calc_munge_mem_gb = calc_munge_memory.rounded_value
    }

    call gsem.gsem_munge as munge {
      input:
        sumstats_files = [preprocess.processed_sumstats[idx]],
        trait_names = [trait_names[idx]],
        sample_sizes = [sample_sizes[idx]],
        ref_snp_list = ref_snp_list,
        out_dir = munge_out_dir,
        info_filter = info_filter,
        maf_filter = maf_filter,
        parallel = parallel,
        cores = cores,
        no_overwrite = no_overwrite,
        ecr_repo = ecr_repo,
        image_source = image_source,
        cpu = munge_cpu,
        mem_gb = select_first([munge_mem_gb, calc_munge_mem_gb])
    }
  }

  if (!defined(ldsc_mem_gb)) {
    call utils.get_total_file_size as ldsc_input_size {
      input:
        input_files = flatten(munge.munged_sumstats)
    }

    Array[String] ldsc_input_dirs = if (ld_dir == wld_dir) then [ld_dir] else [ld_dir, wld_dir]
    call utils.get_total_directory_size as ldsc_dir_size {
      input:
        input_dirs = ldsc_input_dirs
    }

    call utils.round_power_of_two_minus_one as calc_ldsc_memory {
      input:
        input_value = (ldsc_input_size.total_file_size_gb + ldsc_dir_size.total_directory_size_gb) * 1.5
    }
    Int calc_ldsc_mem_gb = calc_ldsc_memory.rounded_value
  }

  call gsem.gsem_ldsc as ldsc {
    input:
      sumstats_files = flatten(munge.munged_sumstats),
      trait_names = trait_names,
      sample_prevs = sample_prevs,
      population_prevs = population_prevs,
      ld_dir = ld_dir,
      wld_dir = wld_dir,
      output_prefix = ldsc_output_prefix,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = ldsc_cpu,
      mem_gb = select_first([ldsc_mem_gb, calc_ldsc_mem_gb])
  }

  if (!defined(commonfactor_mem_gb)) {
    call utils.get_total_file_size as commonfactor_input_size {
      input:
        input_files = [ldsc.ldsc_rds]
    }

    call utils.round_power_of_two_minus_one as calc_commonfactor_memory {
      input:
        input_value = commonfactor_input_size.total_file_size_gb * 1.5
    }
    Int calc_commonfactor_mem_gb = calc_commonfactor_memory.rounded_value
  }

  call gsem.gsem_commonfactor as commonfactor {
    input:
      ldsc_rds = ldsc.ldsc_rds,
      estimation_method = estimation_method,
      output_prefix = commonfactor_output_prefix,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = commonfactor_cpu,
      mem_gb = select_first([commonfactor_mem_gb, calc_commonfactor_mem_gb])
  }

  output {
    Array[File] preprocessed_sumstats = preprocess.processed_sumstats
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File ldsc_rds = ldsc.ldsc_rds
    File ldsc_log = ldsc.ldsc_log
    File commonfactor_rds = commonfactor.commonfactor_rds
  }
}
