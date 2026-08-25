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

workflow gsem_user_model_gwas {
  input {
    Array[File] sumstats_files
    Array[String] trait_names
    Array[Int] sample_sizes
    File ref_snp_list

    Array[Float] sample_prevs
    Array[Float] population_prevs
    String ld_dir
    String wld_dir

    File model_lavaan
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

    Boolean not_printwarn = false
    Array[String]? sub
    Float? toler
    Float? snpse
    String gc = "standard"
    Boolean mpi = false
    Boolean smooth_check = false
    Boolean twas = false
    Boolean std_lv = false
    Boolean not_fix_measurement = false
    Boolean q_snp = false

    String munge_out_dir = "munge_out"
    String ldsc_output_prefix = "ldsc/gsem_ldsc_output"
    String sumstats_output_prefix = "sumstats/gsem_sumstats_output"
    String usergwas_output_prefix = "usergwas/gsem_usergwas_output"

    String? ecr_repo
    String image_source = "docker"
    Int munge_cpu = 1
    String? munge_memory
    Int ldsc_cpu = 1
    String ldsc_memory = "8G"
    Int sumstats_cpu = 1
    String sumstats_memory = "8G"
    Int usergwas_cpu = 1
    String usergwas_memory = "8G"
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
    model_lavaan: "Lavaan model file for userGWAS."
    estimation_method: "Estimation method for userGWAS (for example, DWLS)."
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
    not_printwarn: "If true, suppress per-SNP lavaan warnings/errors in userGWAS output."
    sub: "Optional subset of model components to return from userGWAS."
    toler: "Optional matrix inversion tolerance for userGWAS."
    snpse: "Optional SNP standard error override for userGWAS."
    gc: "Genomic control mode: standard, conserv, or none."
    mpi: "Enable MPI/multi-node mode for userGWAS."
    smooth_check: "Enable smoothing diagnostics in userGWAS."
    twas: "Enable TWAS mode in userGWAS."
    std_lv: "If true, set latent variables to unit variance in userGWAS."
    not_fix_measurement: "If true, do not fix measurement model across SNPs."
    q_snp: "If true, compute Q_SNP statistics."
    munge_out_dir: "Output directory for munged summary statistics."
    ldsc_output_prefix: "Output prefix for LDSC outputs."
    sumstats_output_prefix: "Output prefix for sumstats output."
    usergwas_output_prefix: "Output prefix for user-model GWAS output."
    ecr_repo: "Optional ECR repository URI prefix used with the task ecr image name."
    image_source: "Container source selector: docker or ecr."
    munge_cpu: "Requested CPU cores for munge task runtime."
    munge_memory: "Optional override for munge task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size)."
    ldsc_cpu: "Requested CPU cores for ldsc task runtime."
    ldsc_memory: "Requested memory for ldsc task runtime."
    sumstats_cpu: "Requested CPU cores for sumstats task runtime."
    sumstats_memory: "Requested memory for sumstats task runtime."
    usergwas_cpu: "Requested CPU cores for usergwas task runtime."
    usergwas_memory: "Requested memory for usergwas task runtime."
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

  call gsem.gsem_usergwas as usergwas {
    input:
      ldsc_rds = ldsc.ldsc_rds,
      sumstats_rds = sumstats.sumstats_rds,
      model_lavaan = model_lavaan,
      estimation_method = estimation_method,
      output_prefix = usergwas_output_prefix,
      not_printwarn = not_printwarn,
      sub = sub,
      toler = toler,
      snpse = snpse,
      gc = gc,
      mpi = mpi,
      smooth_check = smooth_check,
      twas = twas,
      std_lv = std_lv,
      not_fix_measurement = not_fix_measurement,
      q_snp = q_snp,
      parallel = parallel,
      cores = cores,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = usergwas_cpu,
      memory = usergwas_memory
  }

  output {
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File ldsc_rds = ldsc.ldsc_rds
    File ldsc_log = ldsc.ldsc_log
    File sumstats_rds = sumstats.sumstats_rds
    File usergwas_rds = usergwas.usergwas_rds
  }
}
