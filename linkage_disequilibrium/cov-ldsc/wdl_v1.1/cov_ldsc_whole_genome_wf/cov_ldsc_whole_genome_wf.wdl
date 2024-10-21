version 1.1

import "cov_ldsc_chr_wf.wdl" as COV_LDSC
import "plink.wdl" as PLINK
import "utils.wdl" as UTILS

workflow cov_ldsc_whole_genome_wf{

    input {

        File bed_in
        File bim_in
        File fam_in
        File cov_eigenvec
        String output_basename

         # Container
        String image_source = "docker"
        String? ecr_repo
       
    }

    call get_chrs {
        input:
            bim_in = bim_in,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get sample count
    call UTILS.wc as sample_count {
        input:
            input_file = fam_in,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int sample_mem_multiplier = floor((sample_count.num_lines / 500) + 1)

    # Get variant count
    call UTILS.wc as variant_count {
        input:
            input_file = bim_in,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int variant_mem_multiplier = floor((variant_count.num_lines / 5000000) + 1)

    scatter(i in range(length(get_chrs.chrs))) {

        String chr = "~{get_chrs.chrs[i]}"

        call PLINK.make_bed as split_dataset {
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                chr = chr,
                output_basename = "~{output_basename}_chr~{chr}",
                cpu = 1,
                mem_gb = sample_mem_multiplier * variant_mem_multiplier,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        call COV_LDSC.cov_ldsc_chr_wf as cov_ldsc {
            input:
                bed_in = split_dataset.bed_out,
                bim_in = split_dataset.bim_out,
                fam_in = split_dataset.fam_out,
                cov_eigenvec = cov_eigenvec,
                output_basename = "~{output_basename}_chr~{chr}",
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    output {
        Array[File] m_files = cov_ldsc.m_file
        Array[File] m_5_50_files = cov_ldsc.m_5_50_file
        Array[File] ldscores = cov_ldsc.ldscore_out
        Array[File] logs = cov_ldsc.log
    }

}

task get_chrs{

    input {
        File bim_in

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 2
    }

    command {
        set -e
        perl -lane 'print $F[0];' ~{bim_in} | sort -u > "chrs.txt"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        Array[String] chrs = read_lines("chrs.txt")
    }

}
