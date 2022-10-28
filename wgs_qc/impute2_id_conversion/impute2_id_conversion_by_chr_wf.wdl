import "biocloud_gwas_workflows/genotype_array_qc/impute2_id_conversion/impute2_id_conversion_wf.wdl" as IDCONVERT

workflow impute2_id_conversion_by_chr_wf{
    Array[File] bed_ins
    Array[File] bim_ins
    Array[File] fam_ins
    Array[String] chrs
    Array[File] id_legend_files
    Array[File] output_basenames

    String build_code

    Int plink_cpu = 1
    Int plink_mem_gb = 2
    Int id_convert_cpu = 1
    Int id_convert_mem_gb = 3

    scatter(index in range(length(chrs))){

        # Convert variant IDs to impute2 format and remove duplicate variants
        call IDCONVERT.impute2_id_conversion_wf as convert_impute2_ids{
            input:
                bed_in = bed_ins[index],
                bim_in = bim_ins[index],
                fam_in = fam_ins[index],
                output_basename = output_basenames[index],
                chrs = chrs[index],
                id_legend_files = id_legend_files[index],
                in_monomorphic_allele = in_monomorphic_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = id_legend_deletion_allele,
                rescue_rsids = rescue_monomorph_rsids,
                remove_duplicates = true,
                build_code = build_code,
                id_convert_cpu = id_convert_cpu,
                id_convert_mem_gb = id_convert_mem_gb,
                plink_cpu = plink_filter_cpu,
                plink_mem_gb = plink_filter_mem_gb
        }
    }

}
