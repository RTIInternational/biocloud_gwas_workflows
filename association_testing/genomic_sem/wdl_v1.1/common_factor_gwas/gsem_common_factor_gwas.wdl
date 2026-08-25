version 1.1

import "genomic_sem.wdl" as gsem

task gsem_calc_munge_memory {
  input {
    File sumstats_file
    File ref_snp_list
  }

  command <<<
    set -euo pipefail

    sumstats_bytes=$(stat -c%s "~{sumstats_file}")
    ref_bytes=$(stat -c%s "~{ref_snp_list}")

    awk -v s="$sumstats_bytes" -v r="$ref_bytes" 'BEGIN {
      mem_gb = ((s + r) / (1024 * 1024 * 1024)) * 1.5;
      if (mem_gb < 1.0) mem_gb = 1.0;
      printf "%.2fG\n", mem_gb;
    }' > munge_memory.txt
  >>>

  output {
    String memory = read_string("munge_memory.txt")
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "1G"
  }
}

workflow gsem_common_factor_gwas {
  input {
    Array[File] sumstats_files
    Array[String] trait_names
    Array[Int] sample_sizes
    File ref_snp_list

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

    String? ecr_repo
    String image_source = "docker"
    Int munge_cpu = 1
    String? munge_memory
    Int ldsc_cpu = 1
    String ldsc_memory = "8G"
    Int sumstats_cpu = 1
    String sumstats_memory = "8G"
    Int commonfactorgwas_cpu = 1
    String commonfactorgwas_memory = "8G"
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
    ldsc_output_prefix: "Output prefix for LDSC outputs."
    sumstats_output_prefix: "Output prefix for sumstats output."
    commonfactorgwas_output_prefix: "Output prefix for common-factor GWAS output."
    ecr_repo: "Optional ECR repository URI prefix used with the task ecr image name."
    image_source: "Container source selector: docker or ecr."
    munge_cpu: "Requested CPU cores for munge task runtime."
    munge_memory: "Optional override for munge task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size)."
    ldsc_cpu: "Requested CPU cores for ldsc task runtime."
    ldsc_memory: "Requested memory for ldsc task runtime."
    sumstats_cpu: "Requested CPU cores for sumstats task runtime."
    sumstats_memory: "Requested memory for sumstats task runtime."
    commonfactorgwas_cpu: "Requested CPU cores for commonfactorgwas task runtime."
    commonfactorgwas_memory: "Requested memory for commonfactorgwas task runtime."
}

  scatter (idx in range(length(sumstats_files))) {
    if (!defined(munge_memory)) {
      call gsem_calc_munge_memory as calc_munge_memory {
        input:
          sumstats_file = sumstats_files[idx],
          ref_snp_list = ref_snp_list
      }
    }

    call gsem.gsem_munge as munge {
      input:
        sumstats_files = [sumstats_files[idx]],
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
        memory = select_first([munge_memory, calc_munge_memory.memory])
    }
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
      memory = ldsc_memory
  }

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
      memory = sumstats_memory
  }

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
      memory = commonfactorgwas_memory
  }

  output {
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File ldsc_rds = ldsc.ldsc_rds
    File ldsc_log = ldsc.ldsc_log
    File sumstats_rds = sumstats.sumstats_rds
    File commonfactorgwas_rds = commonfactorgwas.commonfactorgwas_rds
  }
}
