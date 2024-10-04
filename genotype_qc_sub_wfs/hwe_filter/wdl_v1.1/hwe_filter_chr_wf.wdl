version 1.1

import "plink.wdl" as PLINK
import "normalize_sex_chr_wf.wdl" as NORM

workflow hwe_filter_chr_wf{

    input {

        File bed_in
        File bim_in
        File fam_in
        String output_basename
        String chr
        String build_code

        # HWE p-value filter threshold
        Float hwe_filter_pvalue
        String? hwe_mode

        # Whether or not to include non-founders (i.e. samples with parent_id set) in HWE calculations (default: False)
        Boolean? nonfounders

        # Optional file of individuals to keep
        File? keep_samples

        # Optional file of related individuals to exclude from HWE calculation
        File? related_samples

        # Runtime
        Int plink_chr_cpu = 1
        Int plink_chr_mem_gb = 2
        Int plink_filter_cpu = 1
        Int plink_filter_mem_gb = 2

        String image_source = "docker"
        String? ecr_repo

    }

    # Split PAR/NONPAR if not already split
    if (chr == "23"){
        call NORM.normalize_sex_chr_wf{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                keep_samples = keep_samples,
                expected_chrs = [chr],
                output_basename = "~{output_basename}.sex_norm",
                build_code = build_code,
                no_fail = true,
                plink_cpu = plink_filter_cpu,
                plink_mem_gb = plink_filter_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    File bed_use = select_first([normalize_sex_chr_wf.bed_out, bed_in])
    File bim_use = select_first([normalize_sex_chr_wf.bim_out, bim_in])
    File fam_use = select_first([normalize_sex_chr_wf.fam_out, fam_in])

    # HWE filter
    call PLINK.hardy as hwe{
        input:
            bed_in = bed_use,
            bim_in = bim_use,
            fam_in = fam_use,
            output_basename = "~{output_basename}.chr.~{chr}",
            chrs = [chr],
            filter_females = (chr == "23"),
            hwe_pvalue = hwe_filter_pvalue,
            hwe_mode = hwe_mode,
            nonfounders = nonfounders,
            keep_samples = keep_samples,
            remove_samples = related_samples,
            cpu = plink_chr_cpu,
            mem_gb = plink_chr_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Remove failed HWE SNPs
    call PLINK.make_bed as filter_hwe{
        input:
            bed_in = bed_use,
            bim_in = bim_use,
            fam_in = fam_use,
            chr = chr,
            exclude = hwe.remove,
            keep_samples = keep_samples,
            output_basename = "~{output_basename}.hwe",
            cpu = plink_filter_cpu,
            mem_gb = plink_filter_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        # File bed_out = bed_use
        # File bim_out = bim_use
        # File fam_out = fam_use
        File bed_out = filter_hwe.bed_out
        File bim_out = filter_hwe.bim_out
        File fam_out = filter_hwe.fam_out
    }

}