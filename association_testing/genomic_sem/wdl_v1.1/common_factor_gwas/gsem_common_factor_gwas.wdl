version 1.1

import "../../../../biocloud_wdl_tools/genomic_sem/wdl_v1.1/genomic_sem.wdl" as gsem
import "../../../../biocloud_wdl_tools/utils/wdl_v1.1/utils.wdl" as utils
import "../../../../biocloud_wdl_tools/genomic_sem_preprocessing/wdl_v1.1/genomic_sem_preprocessing.wdl" as preprocessing

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

workflow gsem_common_factor_gwas {
  input {
    Array[File] sumstats_files
    Array[String] trait_names
    Array[Int] sample_sizes
    File ref_snp_list

    SUMSTATS_COLUMNS sumstats_columns

    Array[Float] sample_prevs
    Array[Float] population_prevs
    String ld_dir
    String wld_dir

    String estimation_method = "DWLS"

    Array[Boolean] se_logit
    Array[String]? ols
    Array[String]? linprob
    Array[String]? betas

    Float info_filter = 0.9
    Float maf_filter = 0.01
    Boolean parallel = false
    Int cores = 1
    Boolean no_overwrite = false
    Boolean keep_indel = false

    Float? toler
    Float? snpse
    String gc = "standard"
    Boolean mpi = false
    Boolean smooth_check = false
    Boolean twas = false

    String munge_out_dir = "munge_out"
    String ldsc_output_prefix = "ldsc/gsem_ldsc_output"
    String sumstats_output_prefix = "sumstats/gsem_sumstats_output"
    String commonfactorgwas_output_prefix = "commonfactorgwas/gsem_commonfactorgwas_output"
    String preprocessing_out_dir = "preprocessing_out"

    String? ecr_repo
    String image_source = "docker"

    Int preprocessing_cpu = 1
    Int preprocessing_mem_gb = 4
    Int munge_cpu = 1
    Int? munge_mem_gb
    Int ldsc_cpu = 1
    Int? ldsc_mem_gb
    Int sumstats_cpu = 1
    Int? sumstats_mem_gb
    Int commonfactorgwas_cpu = 1
    Int? commonfactorgwas_mem_gb
}

  parameter_meta {
    sumstats_files: "Input GWAS summary statistics files."
    trait_names: "Trait labels corresponding to sumstats_files."
    sample_sizes: "Per-trait sample sizes (N) for munging/sumstats."
    ref_snp_list: "Reference SNP list (HM3) for munging and sumstats."
    sample_prevs: "Sample prevalences per trait for LDSC."
    population_prevs: "Population prevalences per trait for LDSC."
    ld_dir: "Directory containing LD scores."
    wld_dir: "Directory containing LD score weights."
    estimation_method: "Estimation method for common-factor GWAS (for example, DWLS)."
    se_logit: "Per-trait indicator that SE is on the logistic scale."
    ols: "Optional per-trait OLS indicator for sumstats()."
    linprob: "Optional per-trait linear probability indicator for sumstats()."
    betas: "Optional per-trait beta column names for sumstats()."
    info_filter: "INFO/R2 threshold used by gsem_munge and gsem_sumstats."
    maf_filter: "MAF threshold used by gsem_munge and gsem_sumstats."
    parallel: "Enable parallel processing where supported."
    cores: "Number of cores when parallel processing is enabled."
    no_overwrite: "If true, do not overwrite existing munged outputs."
    keep_indel: "If true, keep indels in gsem_sumstats processing."
    toler: "Optional matrix inversion tolerance for commonfactorgwas."
    snpse: "Optional SNP standard error override for commonfactorgwas."
    gc: "Genomic control mode: standard, conserv, or none."
    mpi: "Enable MPI/multi-node mode for commonfactorgwas."
    smooth_check: "Enable smoothing diagnostics in commonfactorgwas."
    twas: "Enable TWAS mode in commonfactorgwas."
    munge_out_dir: "Output directory for munged summary statistics."
    preprocessing_out_dir: "Output directory for preprocessed summary statistics."
    ldsc_output_prefix: "Output prefix for LDSC outputs."
    sumstats_output_prefix: "Output prefix for sumstats output."
    commonfactorgwas_output_prefix: "Output prefix for common-factor GWAS output."
    ecr_repo: "Optional ECR repository URI prefix used with the task ecr image name."
    image_source: "Container source selector: docker or ecr."
    preprocessing_cpu: "Requested CPU cores for preprocessing task runtime."
    preprocessing_mem_gb: "Requested memory for preprocessing task runtime in GB."
    munge_cpu: "Requested CPU cores for munge task runtime."
    munge_mem_gb: "Optional override for munge task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size), rounded to (2^n - 1) GB."
    ldsc_cpu: "Requested CPU cores for ldsc task runtime."
    ldsc_mem_gb: "Optional override for ldsc task runtime memory; if unset, memory is computed as 1.5x(munged sumstats + ld_dir + wld_dir sizes), rounded to (2^n - 1) GB."
    sumstats_cpu: "Requested CPU cores for sumstats task runtime."
    sumstats_mem_gb: "Optional override for sumstats task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size), rounded to (2^n - 1) GB."
    commonfactorgwas_cpu: "Requested CPU cores for commonfactorgwas task runtime."
    commonfactorgwas_mem_gb: "Optional override for commonfactorgwas task runtime memory; if unset, memory is computed as 1.5x(ldsc_rds + sumstats_rds sizes), rounded to (2^n - 1) GB."
  }

  # Preprocess summary statistics files in parallel
  scatter (idx in range(length(sumstats_files))) {
    call preprocessing.genomic_sem_preprocessing as preprocess {
      input:
        sumstats_file = sumstats_files[idx],
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
        out_file = preprocessing_out_dir + "/preprocessed_" + idx,
        docker_image = "genomic-sem-preprocessing:v1",
        ecr_image = "genomic-sem-preprocessing:v1",
        ecr_repo = ecr_repo,
        image_source = image_source,
        cpu = preprocessing_cpu,
        mem_gb = preprocessing_mem_gb
    }
  }

  scatter (idx in range(length(sumstats_files))) {
    # Calculate the memory requirement for the munge task if it is not defined.
    if (!defined(munge_mem_gb)) {
      call utils.get_total_file_size as munge_input_size {
        input:
          input_files = [preprocess[idx].processed_sumstats, ref_snp_list]
      }
      call utils.round_power_of_two_minus_one as calc_munge_mem {
        input:
          input_value = munge_input_size.total_file_size_gb * 1.5
      }
      Int calc_munge_mem_gb = calc_munge_mem.rounded_value
    } else {
      Int calc_munge_mem_gb = 0
    }

    # Run the gsem_munge task with the calculated or provided memory.
    call gsem.gsem_munge as munge {
      input:
        sumstats_files = [preprocess[idx].processed_sumstats],
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

  # Calculate the memory requirement for the ldsc task if it is not defined.
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
    call utils.round_power_of_two_minus_one as calc_ldsc_mem {
      input:
        input_value = (ldsc_input_size.total_file_size_gb + ldsc_dir_size.total_directory_size_gb) * 1.5
    }
    Int calc_ldsc_mem_gb = calc_ldsc_mem.rounded_value
  } else {
    Int calc_ldsc_mem_gb = 0
  }

  # Run the ldsc task with the calculated or provided memory.
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

  # Calculate the memory requirement for the sumstats task if it is not defined.
  if (!defined(sumstats_mem_gb)) {
    call utils.get_total_file_size as sumstats_input_size {
      input:
        input_files = flatten([sumstats_files, ref_snp_list])
    }
    call utils.round_power_of_two_minus_one as calc_sumstats_mem {
      input:
        input_value = sumstats_input_size.total_file_size_gb * 1.5
    }
    Int calc_sumstats_mem_gb = calc_sumstats_mem.rounded_value
  } else {
    Int calc_sumstats_mem_gb = 0
  }

  # Run the sumstats task with the calculated or provided memory.
  call gsem.gsem_sumstats as sumstats {
    input:
      sumstats_files = sumstats_files,
      trait_names = trait_names,
      sample_sizes = sample_sizes,
      se_logit = se_logit,
      ref_snp_list = ref_snp_list,
      output_prefix = sumstats_output_prefix,
      ols = ols,
      linprob = linprob,
      betas = betas,
      info_filter = info_filter,
      maf_filter = maf_filter,
      keep_indel = keep_indel,
      parallel = parallel,
      cores = cores,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = sumstats_cpu,
      mem_gb = select_first([sumstats_mem_gb, calc_sumstats_mem_gb])
  }

  # Calculate the memory requirement for the commonfactorgwas task if it is not defined.
  if (!defined(commonfactorgwas_mem_gb)) {
    call utils.get_total_file_size as commonfactorgwas_input_size {
      input:
        input_files = [ldsc.ldsc_rds, sumstats.sumstats_rds]
    }
    call utils.round_power_of_two_minus_one as calc_commonfactorgwas_mem {
      input:
        input_value = commonfactorgwas_input_size.total_file_size_gb * 1.5
    }
    Int calc_commonfactorgwas_mem_gb = calc_commonfactorgwas_mem.rounded_value
  } else {
    Int calc_commonfactorgwas_mem_gb = 0
  }

  # Run the commonfactorgwas task with the calculated or provided memory.
  call gsem.gsem_commonfactorgwas as commonfactorgwas {
    input:
      ldsc_rds = ldsc.ldsc_rds,
      sumstats_rds = sumstats.sumstats_rds,
      estimation_method = estimation_method,
      output_prefix = commonfactorgwas_output_prefix,
      toler = toler,
      snpse = snpse,
      gc = gc,
      mpi = mpi,
      smooth_check = smooth_check,
      twas = twas,
      parallel = parallel,
      cores = cores,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = commonfactorgwas_cpu,
      mem_gb = select_first([commonfactorgwas_mem_gb, calc_commonfactorgwas_mem_gb])
  }

  output {
    Array[File] preprocessed_sumstats = preprocess.processed_sumstats
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File ldsc_rds = ldsc.ldsc_rds
    File ldsc_log = ldsc.ldsc_log
    File sumstats_rds = sumstats.sumstats_rds
    File commonfactorgwas_rds = commonfactorgwas.commonfactorgwas_rds
  }
}
