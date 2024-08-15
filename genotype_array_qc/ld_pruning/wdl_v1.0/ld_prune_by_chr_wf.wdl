version 1.0

import "ld_prune_wf.wdl" as LD

workflow ld_prune_by_chr_wf{

    input {

        Array[File] bed_ins
        Array[File] bim_ins
        Array[File] fam_ins
        Array[String] output_basenames
        Array[String] chrs

        File? exclude_regions
        Float? maf
        String ld_type = "indep-pairphase"
        Int? window_size
        Int? step_size
        String? window_size_unit
        Float? r2_threshold
        Float? vif_threshold
        Int? x_chr_mode

        Int cpu
        Int mem_gb

        # Container
        String image_source = "docker"
        String? ecr_repo

    }

    scatter(index in range(length(chrs))){

        call LD.ld_prune_wf as ld_prune{
            input:
                bed_in = bed_ins[index],
                bim_in = bim_ins[index],
                fam_in = fam_ins[index],
                output_basename = output_basenames[index],
                ld_type = ld_type,
                window_size = window_size,
                step_size = step_size,
                r2_threshold = r2_threshold,
                cpu = cpu,
                mem_gb = mem_gb,
                maf = maf,
                chr = chrs[index],
                exclude_regions = exclude_regions,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

}
