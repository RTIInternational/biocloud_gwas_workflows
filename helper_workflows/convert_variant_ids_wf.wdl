import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT

workflow convert_variant_ids_wf{
    Array[File] input_files
    Array[File] ref_files
    Array[String] chrs
    Array[File] output_files

    Int in_header = 0
    String in_sep
    Int in_id_col
    Int in_chr_col
    Int in_pos_col
    Int in_a1_col
    Int in_a2_col
    String in_missing_allele
    String in_deletion_allele
    String ref_deletion_allele
    Int chunk_size
    Boolean? rescue_rsids
    String? output_compression

    # Resources
    Int cpu = 1
    Int mem_gb = 3


    # Parallelize
    scatter(i in range(length(input_files))){

        # Convert IDs
        call IDCONVERT.convert_variant_ids{
            input:
                in_file = input_files[i],
                ref = ref_files[i],
                chr = chrs[i],
                in_header = in_header,
                in_sep = in_sep,
                in_id_col = in_id_col,
                in_chr_col = in_chr_col,
                in_pos_col = in_pos_col,
                in_a1_col = in_a1_col,
                in_a2_col = in_a2_col,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = ref_deletion_allele,
                chunk_size = chunk_size,
                output_filename = output_files[i],
                rescue_rsids = if defined(rescue_rsids) then rescue_rsids else null,
                output_compression = if defined(output_compression) then output_compression else null,
                cpu = cpu,
                mem_gb = mem_gb
        }

    }

    output{
        Array[File] output_files = output_files
    }
}