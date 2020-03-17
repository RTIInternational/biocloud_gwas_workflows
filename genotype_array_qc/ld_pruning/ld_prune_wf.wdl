import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow ld_prune_wf{
    File bed_in
    File bim_in
    File fam_in
    File output_basename

    File? exclude_regions
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
            exclude_regions = exclude_regions
    }

    # Filter to include only the LD-pruned markers returned from previous step
    call PLINK.make_bed{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            output_basename = output_basename,
            chr = chr,
            extract = prune_ld_markers.include_markers
    }

    output{
        File bed_out = make_bed.bed_out
        File bim_out = make_bed.bim_out
        File fam_out = make_bed.fam_out
        File include_markers = prune_ld_markers.include_markers
        File exclude_markers = prune_ld_markers.exclude_markers
    }
}