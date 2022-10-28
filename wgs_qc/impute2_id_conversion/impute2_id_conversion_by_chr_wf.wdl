import "biocloud_gwas_workflows/genotype_array_qc/impute2_id_conversion/impute2_id_conversion_wf.wdl" as IMPUTE2

workflow impute2_id_conversion_by_chr_wf{
    Array[File] bed_ins
    Array[File] bim_ins
    Array[File] fam_ins
    Array[String] chrs
    Array[File] id_legend_files
    Array[File] output_basenames

    String build_code

    Int plink_cpu
    Int plink_mem_gb
    Int id_convert_cpu
    Int id_convert_mem_gb

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
                exclude_regions = exclude_regions
        }
    }

}
