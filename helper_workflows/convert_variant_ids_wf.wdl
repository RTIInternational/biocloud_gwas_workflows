import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT
import "biocloud_gwas_workflows/biocloud_wdl_tools/rti-tsv-utils/rti-tsv-utils.wdl" as SORT

workflow convert_variant_ids_wf{
    Array[File] input_files
    Array[File] ref_files
    Array[String] chrs
    Array[File] output_files

    Int in_header
    String in_sep
    Int in_id_col
    Int in_chr_col
    Int in_pos_col
    Int in_a1_col
    Int in_a2_col
    String in_missing_allele
    String in_deletion_allele
    String ref_deletion_allele
    Int? in_chunk_size
    Int? ref_chunk_size
    Boolean? rescue_rsids
    String? output_compression

    Boolean sort = false
    String sort_col = "pos"

    # Convert resources
    Int cpu = 1
    Int mem_gb = 4
    String docker = "rtibiocloud/convert_variant_ids:v1_7b4adf9"

    # Sort resources
    Int sort_cpu = 1
    Int sort_mem_gb = 4
    String? sort_docker
   
    # Parallelize
    scatter(i in range(length(input_files))){

        # Sort
        if(sort){
            String new_sep = "tab"
            call SORT.tsv_sort as tsv_sort{
                input:
                    in_file = input_files[i],
                    cols = sort_col,
                    out_prefix = "sort_chr${chrs[i]}",
                    in_file_sep = in_sep,
                    cpu = sort_cpu,
                    mem_gb = sort_mem_gb,
                    docker = sort_docker
            }
        }

        File in_file = select_first([tsv_sort.out_tsv, input_files[i]])
        String sep = select_first([new_sep, in_sep])

        # Convert IDs
        call IDCONVERT.convert_variant_ids{
            input:
                chr = chrs[i],
                in_file = in_file,
                in_header = in_header,
                in_sep = sep,
                in_id_col = in_id_col,
                in_chr_col = in_chr_col,
                in_pos_col = in_pos_col,
                in_a1_col = in_a1_col,
                in_a2_col = in_a2_col,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                in_chunk_size = in_chunk_size,
                ref = ref_files[i],
                ref_deletion_allele = ref_deletion_allele,
                ref_chunk_size = ref_chunk_size,
                output_filename = output_files[i],
                rescue_rsids = rescue_rsids,
                output_compression = output_compression,
                cpu = cpu,
                mem_gb = mem_gb,
                docker = docker
        }

    }

    output{
        Array[File] results_files = output_files
    }
}
