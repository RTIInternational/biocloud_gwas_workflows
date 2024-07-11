version 1.1

import "plink.wdl" as PLINK
import "utils.wdl" as UTILS

workflow ld_prune_wf{

    input {

        File bed_in
        File bim_in
        File fam_in
        File output_basename

        File? exclude_regions
        File? exclude
        String? chr
        Float? maf
        String ld_type = "indep-pairphase"
        Int window_size
        Int step_size
        String? window_size_unit
        Float? r2_threshold
        Float? vif_threshold
        Int? x_chr_mode

        Int cpu
        Int mem_gb

        # Container
        String container_source = "docker"
        Int? ecr_account_id

    }

    # Count number of samples (if <50 set --bad-ld option or plink2 will error out
    call UTILS.wc as get_num_samples{
        input:
            input_file = fam_in,
            container_source = container_source,
            ecr_account_id = ecr_account_id
    }

    # Get LD pruning set
    call PLINK.prune_ld_markers{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            output_basename = output_basename,
            ld_type = ld_type,
            window_size = window_size,
            step_size = step_size,
            window_size_unit = window_size_unit,
            r2_threshold = r2_threshold,
            vif_threshold = vif_threshold,
            x_chr_mode = x_chr_mode,
            cpu = cpu,
            mem_gb = mem_gb,
            maf = maf,
            chr = chr,
            bad_ld = (get_num_samples.num_lines <= 50),
            exclude_regions = exclude_regions,
            exclude = exclude,
            container_source = container_source,
            ecr_account_id = ecr_account_id
    }

    # Filter to include only the LD-pruned markers returned from previous step
    call PLINK.make_bed{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            output_basename = output_basename,
            chr = chr,
            cpu = cpu,
            mem_gb = mem_gb,
            extract = prune_ld_markers.include_markers,
            container_source = container_source,
            ecr_account_id = ecr_account_id
    }

    output{
        File bed_out = make_bed.bed_out
        File bim_out = make_bed.bim_out
        File fam_out = make_bed.fam_out
        File include_markers = prune_ld_markers.include_markers
        File exclude_markers = prune_ld_markers.exclude_markers
    }
}
