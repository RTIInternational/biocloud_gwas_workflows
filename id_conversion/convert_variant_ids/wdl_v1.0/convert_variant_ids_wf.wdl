version 1.0

import "convert_variant_ids.wdl" as IDCONVERT

workflow convert_variant_ids_wf{

    input {

        # Require parameters
        Array[File] ref_files
        Array[String] chrs
        Array[File] input_files
        Array[File] output_files
        Int in_id_col
        Int in_chr_col
        Int in_pos_col
        Int in_a1_col
        Int in_a2_col

        # Optional parameters
        Int? in_header
        String? in_sep
        String? in_missing_allele
        String? in_deletion_allele
        String? ref_deletion_allele
        Int? in_chunk_size
        Int? ref_chunk_size
        Boolean? rescue_rsids
        String? output_compression

        # Convert resources
        Int convert_cpu = 2
        Int convert_mem_gb = 4

        String image_source = "docker"
        String? ecr_repo

    }

    scatter(i in range(length(chrs))) {

        # Convert IDs
        call IDCONVERT.convert_variant_ids{
            input:
                ref = ref_files[i],
                chr = chrs[i],
                in_file = input_files[i],
                output_filename = output_files[i],
                in_id_col = in_id_col,
                in_chr_col = in_chr_col,
                in_pos_col = in_pos_col,
                in_a1_col = in_a1_col,
                in_a2_col = in_a2_col,
                in_header = in_header,
                in_sep = in_sep,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                in_chunk_size = in_chunk_size,
                ref_deletion_allele = ref_deletion_allele,
                ref_chunk_size = ref_chunk_size,
                rescue_rsids = rescue_rsids,
                output_compression = output_compression,
                cpu = convert_cpu,
                mem_gb = convert_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    output{
        Array[File] results_files = convert_variant_ids.output_file
        Array[File] log_files = convert_variant_ids.log
    }
}
