import "biocloud_gwas_workflows/genotype_array_qc/impute2_id_conversion/impute2_id_conversion_wf.wdl" as IDCONVERT

workflow impute2_id_conversion_by_chr_wf{
    Array[File] bed_ins
    Array[File] bim_ins
    Array[File] fam_ins
    Array[String] chrs
    Array[File] id_legend_files
    Array[File] output_basenames

    String build_code
    String in_monomorphic_allele = "0"
    String in_deletion_allele = "-"
    String id_legend_deletion_allele = "."
    Boolean rescue_monomorph_rsids = true

    Int plink_cpu = 1
    Int plink_mem_gb = 2
    Int id_convert_cpu = 1
    Int id_convert_mem_gb = 3

    scatter(i in range(length(chrs))){
        Int chr = chrs[i]

        Array[String] chr_array = [chr]
        # Convert variant IDs to impute2 format and remove duplicate variants
        call IDCONVERT.impute2_id_conversion_wf as convert_impute2_ids{
            input:
                bed_in = bed_ins[i],
                bim_in = bim_ins[i],
                fam_in = fam_ins[i],
                output_basename = output_basenames[i],
                chrs = chr_array,
                id_legend_files = id_legend_files,
                in_monomorphic_allele = in_monomorphic_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = id_legend_deletion_allele,
                rescue_rsids = rescue_monomorph_rsids,
                remove_duplicates = true,
                build_code = build_code,
                id_convert_cpu = id_convert_cpu,
                id_convert_mem_gb = id_convert_mem_gb,
                plink_cpu = plink_cpu,
                plink_mem_gb = plink_mem_gb
        }
    }

}
