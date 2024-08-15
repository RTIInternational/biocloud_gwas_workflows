version 1.0

import "plink.wdl" as PLINK
import "convert_variant_ids.wdl" as IDCONVERT
import "rti-tsv-utils.wdl" as SORT

workflow convert_variant_ids_bfile_chr_wf{

    input{

        File bim_in
        File ref_file
        String chr
        String output_filename

        String? in_sep
        String? in_missing_allele
        String? in_deletion_allele
        String? ref_deletion_allele
        Int? in_chunk_size
        Int? ref_chunk_size
        Boolean? rescue_rsids
        String? output_compression

        # Convert resources
        Int convert_cpu = 1
        Int convert_mem_gb = 4

        String image_source = "docker"
        String? ecr_repo

    }

    # Convert IDs
    call IDCONVERT.convert_variant_ids{
        input:
            chr = chr,
            in_file = bim_in,
            in_header = 0,
            in_sep = in_sep,
            in_id_col = 1,
            in_chr_col = 0,
            in_pos_col = 3,
            in_a1_col = 4,
            in_a2_col = 5,
            in_missing_allele = in_missing_allele,
            in_deletion_allele = in_deletion_allele,
            in_chunk_size = in_chunk_size,
            ref = ref_file,
            ref_deletion_allele = ref_deletion_allele,
            ref_chunk_size = ref_chunk_size,
            output_filename = output_filename,
            rescue_rsids = rescue_rsids,
            output_compression = output_compression,
            cpu = convert_cpu,
            mem_gb = convert_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File bim_out = convert_variant_ids.output_file
        File log_out = convert_variant_ids.log
    }
}
