import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD

workflow ld_prune_by_chr_wf{
    Array[File] bed_ins
    Array[File] bim_ins
    Array[File] fam_ins
    Array[File] output_basenames
    Array[String] chrs

    File? exclude_regions
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
                cpu = ld_cpu,
                mem_gb = ld_mem_gb,
                maf = min_ld_maf,
                chr = chrs[index],
                exclude_regions = ld_exclude_regions
        }
    }

}
