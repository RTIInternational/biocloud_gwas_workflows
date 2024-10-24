version 1.1

import "cov_ldsc.wdl" as COV_LDSC
import "utils.wdl" as UTILS

workflow cov_ldsc_chr_wf{

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

    # Get sample count
    call UTILS.wc as sample_count {
        input:
            input_file = fam_in,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get variant count
    call UTILS.wc as variant_count {
        input:
            input_file = bim_in,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    Int mem_gb = floor((sample_count.num_lines * variant_count.num_lines) / 60000000)
    Int cpu = floor(mem_gb / 8)

    call COV_LDSC.cov_ldsc {
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            cov_eigenvec = cov_eigenvec,
            output_basename = output_basename,
            cpu = cpu,
            mem_gb = mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output {
        File m_file = cov_ldsc.m_file
        File m_5_50_file = cov_ldsc.m_5_50_file
        File ldscore_out = cov_ldsc.ldscore_out
        File log = cov_ldsc.log
    }

}