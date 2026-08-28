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

workflow gsem_common_factor_gwas {
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

    Int tsv_append_cpu = 1
    Int tsv_append_mem_gb = 2
    Int gwas_chunk_size = 500000
    Int tsv_split_cpu = 1
    Int tsv_split_mem_gb = 2
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
    Int merge_rds_cpu = 1
    Int merge_rds_mem_gb = 8
}

  parameter_meta {
    sumstats_files: "Input GWAS summary statistics files (either one genome sumstats file or per-chromosome sumstats files per trait/study)."
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
    tsv_append_cpu: "Requested CPU cores for tsv_append task runtime."
    tsv_append_mem_gb: "Memory in GB for tsv_append task runtime."
    gwas_chunk_size: "Number of lines (variants) per split chunk for sumstats and GWAS analysis."
    tsv_split_cpu: "Requested CPU cores for tsv_split task runtime."
    tsv_split_mem_gb: "Memory in GB for tsv_split task runtime."
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
    merge_rds_cpu: "Requested CPU cores for merge_rds task runtime."
    merge_rds_mem_gb: "Memory in GB for merge_rds task runtime."
  }

  Int num_traits = length(sumstats_files)
  Boolean input_lengths_match = (
    length(trait_names) == num_traits &&
    length(sample_sizes) == num_traits &&
    length(sample_prevs) == num_traits &&
    length(population_prevs) == num_traits &&
    length(se_logit) == num_traits &&
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
        msg = "Workflow input error: sumstats_files, trait_names, sample_sizes, sample_prevs, population_prevs, se_logit, and all arrays within sumstats_columns must have the same length.",
        ecr_repo = ecr_repo,
        image_source = image_source
    }
  }

  # Preprocess summary statistics files in parallel
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
    # Calculate the memory requirement for the munge task if it is not defined.
    if (!defined(munge_mem_gb)) {
      call utils.get_total_file_size as munge_input_size {
        input:
          input_files = [preprocess.processed_sumstats[idx], ref_snp_list]
      }
      call utils.round_power_of_two_minus_one as calc_munge_mem {
        input:
          input_value = munge_input_size.total_file_size_gb * 1.5
      }
      Int calc_munge_mem_gb = calc_munge_mem.rounded_value
    }

    # Run the gsem_munge task with the calculated or provided memory.
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
        input_files = flatten([preprocess.processed_sumstats, [ref_snp_list]])
    }
    call utils.round_power_of_two_minus_one as calc_sumstats_mem {
      input:
        input_value = sumstats_input_size.total_file_size_gb * 1.5
    }
    Int calc_sumstats_mem_gb = calc_sumstats_mem.rounded_value
  }

  # Run the sumstats task with the calculated or provided memory.
  call gsem.gsem_sumstats as sumstats {
    input:
      sumstats_files = preprocess.processed_sumstats,
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

  call rti_tsv.tsv_split as split_sumstats {
    input:
      input_file = sumstats.sumstats_tsv,
      lines_per_split = gwas_chunk_size,
      output_prefix = "~{preprocessing_out_dir}/split_sumstats",
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = tsv_split_cpu,
      mem_gb = tsv_split_mem_gb
  }

  scatter (chunk_idx in range(length(split_sumstats.split_files))) {
    File current_chunk_sumstats = split_sumstats.split_files[chunk_idx]

    # Calculate the memory requirement for the commonfactorgwas task if it is not defined.
    if (!defined(commonfactorgwas_mem_gb)) {
      call utils.get_total_file_size as commonfactorgwas_chunk_input_size {
        input:
          input_files = [ldsc.ldsc_rds, current_chunk_sumstats]
      }
      call utils.round_power_of_two_minus_one as calc_commonfactorgwas_chunk_mem {
        input:
          input_value = commonfactorgwas_chunk_input_size.total_file_size_gb * 1.5
      }
      Int calc_commonfactorgwas_chunk_mem_gb = calc_commonfactorgwas_chunk_mem.rounded_value
    }

    # Run the commonfactorgwas task with the calculated or provided memory.
    call gsem.gsem_commonfactorgwas as commonfactorgwas_chunk {
      input:
        ldsc_rds = ldsc.ldsc_rds,
        sumstats = current_chunk_sumstats,
        estimation_method = estimation_method,
        output_prefix = "~{commonfactorgwas_output_prefix}_chunk_~{chunk_idx}",
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
        mem_gb = select_first([commonfactorgwas_mem_gb, calc_commonfactorgwas_chunk_mem_gb])
    }
  }

  call gsem.gsem_merge_rds as merge_commonfactorgwas {
    input:
      rds_files = commonfactorgwas_chunk.commonfactorgwas_rds,
      output_prefix = commonfactorgwas_output_prefix,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = merge_rds_cpu,
      mem_gb = merge_rds_mem_gb
  }

  output {
    Array[File] preprocessed_sumstats = preprocess.processed_sumstats
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File ldsc_rds = ldsc.ldsc_rds
    File ldsc_log = ldsc.ldsc_log
    File sumstats_rds = sumstats.sumstats_rds
    File commonfactorgwas_rds = merge_commonfactorgwas.merged_rds
  }
}
