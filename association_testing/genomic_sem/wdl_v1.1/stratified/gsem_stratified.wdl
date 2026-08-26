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

workflow gsem_stratified {
  input {
    Array[File] sumstats_files
    Array[String] trait_names
    Array[Int] sample_sizes
    File ref_snp_list

    Array[Float] sample_prevs
    Array[Float] population_prevs
    String ld_dir
    String wld_dir
    String frq_dir

    File model_lavaan
    File params
    File? fixparam

    Float info_filter = 0.9
    Float maf_filter = 0.01
    Boolean parallel = false
    Int cores = 1
    Boolean no_overwrite = false

    Int n_blocks = 200
    Boolean include_cont = false
    String fix = "regressions"
    Boolean std_lv = false
    Boolean not_rm_flank = false
    Boolean tau = false
    Boolean not_base = false
    Float? toler

    String munge_out_dir = "munge_out"
    String s_ldsc_output_prefix = "s_ldsc/gsem_s_ldsc_output"
    String enrich_output_prefix = "enrich/gsem_enrich_output"

    String? ecr_repo
    String image_source = "docker"
    String preprocessing_out_dir = "preprocessing_out"
    Int preprocessing_cpu = 1
    Int preprocessing_mem_gb = 4
    SUMSTATS_COLUMNS sumstats_columns
    Int munge_cpu = 1
    Int? munge_mem_gb
    Int s_ldsc_cpu = 1
    Int? s_ldsc_mem_gb
    Int enrich_cpu = 1
    Int? enrich_mem_gb
}

  parameter_meta {
    sumstats_files: "Input GWAS summary statistics files."
    trait_names: "Trait labels corresponding to sumstats_files."
    sample_sizes: "Per-trait sample sizes (N) for munging."
    ref_snp_list: "Reference SNP list (HM3) for munging."
    sample_prevs: "Sample prevalences per trait for stratified LDSC."
    population_prevs: "Population prevalences per trait for stratified LDSC."
    ld_dir: "Directory containing partitioned LD scores."
    wld_dir: "Directory containing LD score weights."
    frq_dir: "Directory containing allele frequency files for s_ldsc."
    model_lavaan: "Lavaan model file for enrich()."
    params: "File listing target lavaan parameters for enrichment estimation."
    fixparam: "Optional file listing parameters to fix during enrichment estimation."
    info_filter: "INFO/R2 threshold used by gsem_munge."
    maf_filter: "MAF threshold used by gsem_munge."
    parallel: "Enable parallel processing where supported."
    cores: "Number of cores when parallel processing is enabled."
    no_overwrite: "If true, do not overwrite existing munged outputs."
    n_blocks: "Number of jackknife blocks for s_ldsc standard error estimation."
    include_cont: "If true, include continuous annotations in s_ldsc."
    fix: "Parameter class to fix in enrichment model (regressions, variances, covariances)."
    std_lv: "If true, use unit-variance identification for latent variables."
    not_rm_flank: "If true, keep flanking-window and continuous annotations in enrich output."
    tau: "If true, use tau matrices instead of zero-order matrices in enrich."
    not_base: "If true, exclude baseline model estimates from enrich output."
    toler: "Optional matrix inversion tolerance for enrich."
    munge_out_dir: "Output directory for munged summary statistics."
    s_ldsc_output_prefix: "Output prefix for stratified LDSC outputs."
    enrich_output_prefix: "Output prefix for enrich output."
    ecr_repo: "Optional ECR repository URI prefix used with the task ecr image name."
    image_source: "Container source selector: docker or ecr."
    preprocessing_out_dir: "Output directory for preprocessed summary statistics."
    sumstats_columns: "Column mapping for preprocessing."
    preprocessing_cpu: "Requested CPU cores for preprocessing task runtime."
    preprocessing_mem_gb: "Memory in GB for preprocessing task runtime."
    munge_cpu: "Requested CPU cores for munge task runtime."
    munge_mem_gb: "Optional override for munge task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size), rounded to (2^n - 1) GB."
    s_ldsc_cpu: "Requested CPU cores for s_ldsc task runtime."
    s_ldsc_mem_gb: "Optional override for s_ldsc task runtime memory; if unset, memory is computed as 1.5x(munged sumstats + ld_dir + wld_dir + frq_dir sizes), rounded to (2^n - 1) GB."
    enrich_cpu: "Requested CPU cores for enrich task runtime."
    enrich_mem_gb: "Optional override for enrich task runtime memory; if unset, memory is computed as 1.5x(s_ldsc_rds + model_lavaan + params sizes), rounded to (2^n - 1) GB."
}

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
    if (!defined(munge_mem_gb)) {
      call utils.get_total_file_size as munge_input_size {
        input:
          input_files = [preprocess[idx].processed_sumstats, ref_snp_list]
      }

      call utils.round_power_of_two_minus_one as calc_munge_memory {
        input:
          input_value = munge_input_size.total_file_size_gb * 1.5
      }
      Int calc_munge_mem_gb = calc_munge_memory.rounded_value
    } else {
      Int calc_munge_mem_gb = 0
    }

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

  if (!defined(s_ldsc_mem_gb)) {
    call utils.get_total_file_size as s_ldsc_input_size {
      input:
        input_files = flatten(munge.munged_sumstats)
    }

    Array[String] s_ldsc_input_dirs = if (ld_dir == wld_dir) then
      if (frq_dir == ld_dir) then [ld_dir] else [ld_dir, frq_dir]
    else
      if (frq_dir == ld_dir) then [ld_dir, wld_dir] else
      if (frq_dir == wld_dir) then [ld_dir, wld_dir] else [ld_dir, wld_dir, frq_dir]

    call utils.get_total_directory_size as s_ldsc_dir_size {
      input:
        input_dirs = s_ldsc_input_dirs
    }

    call utils.round_power_of_two_minus_one as calc_s_ldsc_memory {
      input:
        input_value = (s_ldsc_input_size.total_file_size_gb + s_ldsc_dir_size.total_directory_size_gb) * 1.5
    }
    Int calc_s_ldsc_mem_gb = calc_s_ldsc_memory.rounded_value
  } else {
    Int calc_s_ldsc_mem_gb = 0
  }

  call gsem.gsem_s_ldsc as s_ldsc {
    input:
      sumstats_files = flatten(munge.munged_sumstats),
      trait_names = trait_names,
      sample_prevs = sample_prevs,
      population_prevs = population_prevs,
      ld_dir = ld_dir,
      wld_dir = wld_dir,
      frq_dir = frq_dir,
      output_prefix = s_ldsc_output_prefix,
      n_blocks = n_blocks,
      include_cont = include_cont,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = s_ldsc_cpu,
      mem_gb = select_first([s_ldsc_mem_gb, calc_s_ldsc_mem_gb])
  }

  if (!defined(enrich_mem_gb)) {
    call utils.get_total_file_size as enrich_input_size {
      input:
        input_files = [s_ldsc.s_ldsc_rds, model_lavaan, params]
    }

    call utils.round_power_of_two_minus_one as calc_enrich_memory {
      input:
        input_value = enrich_input_size.total_file_size_gb * 1.5
    }
    Int calc_enrich_mem_gb = calc_enrich_memory.rounded_value
  } else {
    Int calc_enrich_mem_gb = 0
  }

  call gsem.gsem_enrich as enrich {
    input:
      s_ldsc_rds = s_ldsc.s_ldsc_rds,
      model_lavaan = model_lavaan,
      params = params,
      output_prefix = enrich_output_prefix,
      fix = fix,
      std_lv = std_lv,
      not_rm_flank = not_rm_flank,
      tau = tau,
      not_base = not_base,
      toler = toler,
      fixparam = fixparam,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = enrich_cpu,
      mem_gb = select_first([enrich_mem_gb, calc_enrich_mem_gb])

  output {
    Array[File] preprocessed_sumstats = preprocess.processed_sumstats
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File s_ldsc_rds = s_ldsc.s_ldsc_rds
    File s_ldsc_log = s_ldsc.s_ldsc_log
    File enrich_rds = enrich.enrich_rds
  }
}
