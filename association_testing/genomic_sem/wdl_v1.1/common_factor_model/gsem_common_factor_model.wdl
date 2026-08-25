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

workflow gsem_common_factor_model {
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
    Int munge_cpu = 1
    String? munge_memory
    Int ldsc_cpu = 1
    String ldsc_memory = "8G"
    Int commonfactor_cpu = 1
    String commonfactor_memory = "8G"
}

  parameter_meta {
    sumstats_files: "Input GWAS summary statistics files."
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
    munge_cpu: "Requested CPU cores for munge task runtime."
    munge_memory: "Optional override for munge task runtime memory; if unset, memory is computed as 1.5x(sumstats file size + ref_snp_list size)."
    ldsc_cpu: "Requested CPU cores for ldsc task runtime."
    ldsc_memory: "Requested memory for ldsc task runtime."
    commonfactor_cpu: "Requested CPU cores for commonfactor task runtime."
    commonfactor_memory: "Requested memory for commonfactor task runtime."
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

  call gsem.gsem_commonfactor as commonfactor {
    input:
      ldsc_rds = ldsc.ldsc_rds,
      estimation_method = estimation_method,
      output_prefix = commonfactor_output_prefix,
      ecr_repo = ecr_repo,
      image_source = image_source,
      cpu = commonfactor_cpu,
      memory = commonfactor_memory
  }

  output {
    Array[File] munged_sumstats = flatten(munge.munged_sumstats)
    File ldsc_rds = ldsc.ldsc_rds
    File ldsc_log = ldsc.ldsc_log
    File commonfactor_rds = commonfactor.commonfactor_rds
  }
}
