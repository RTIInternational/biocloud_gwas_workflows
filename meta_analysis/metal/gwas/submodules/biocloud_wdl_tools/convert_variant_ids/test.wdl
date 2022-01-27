import "biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as ID

workflow test_convert_variant_ids{
    File bim
    Array[String] chrs
    Array[File] refs
    String output_basename = "derp"
    Int cpu
    Int mem_gb

    scatter(index in range(length(chrs))){
        String chr = chrs[index]
        call ID.convert_variant_ids{
            input:
                in_file = bim,
                ref = refs[index],
                chr = chrs[index],
                in_header = 0,
                in_sep = "tab",
                in_id_col = 1,
                in_chr_col = 0,
                in_pos_col = 3,
                in_a1_col = 4,
                in_a2_col = 5,
                in_missing_allele = 0,
                in_deletion_allele = "-",
                ref_deletion_allele = ".",
                output_filename = "${output_basename}.${chr}.bim",
                cpu = cpu,
                mem_gb = mem_gb,
                rescue_rsids = true
        }
    }

    output{
        Array[File] converted_id_bims = convert_variant_ids.output_file
        Array[File] converted_id_logs = convert_variant_ids.log
    }
}

