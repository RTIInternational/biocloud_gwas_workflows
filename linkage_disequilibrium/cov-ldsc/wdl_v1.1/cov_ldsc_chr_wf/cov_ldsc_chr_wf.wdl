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
    Int sample_mem_multiplier = floor((sample_count.num_lines / 500) + 1)

    # Get variant count
    call UTILS.wc as variant_count {
        input:
            input_file = bim_in,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int variant_mem_multiplier = floor((variant_count.num_lines / 5000000) + 1)

    call COV_LDSC.cov_ldsc {
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            cov_eigenvec = cov_eigenvec,
            output_basename = output_basename,
            cpu = 2,
            mem_gb = 4 * sample_mem_multiplier * variant_mem_multiplier,
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