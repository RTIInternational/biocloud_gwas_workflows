import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/king/king.wdl" as KING
import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

workflow king_kinship_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    # Degree of relatedness to consider
    Int degree

    # Number of subsets to make for splitting
    Int num_splits = 4

    Int plink_cpu = 1
    Int plink_mem_gb = 2

    Int king_split_cpu = 2
    Int king_split_mem_gb = 4

    # Split plink sample file into however many chunks
    call SPLIT.split_file as split_fam{
        input:
            input_file = fam_in,
            output_basename = basename(fam_in, ".fam"),
            output_extension = ".fam",
            num_splits = num_splits
    }

    # Create split plink files
    scatter(split_index in range(length(split_fam.output_files))){

        # Generate keep file from split fam
        call UTILS.cut{
            input:
                input_file = split_fam.output_files[split_index],
                output_filename = "${output_basename}.split.${split_index}.keep",
                args = "-f 1,2"
        }

        # Subset plink dataset using keep file
        call PLINK.make_bed as split_plink{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                keep_samples = cut.output_file,
                output_basename = "${output_basename}.split.${split_index}",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Do kinship within each subset
        call KING.kinship as subset_kinships{
            input:
                bed_in = split_plink.bed_out,
                bim_in = split_plink.bim_out,
                fam_in = split_plink.fam_out,
                output_basename = "${output_basename}.split.${split_index}",
                degree = degree,
                cpu = king_split_cpu,
                mem_gb = king_split_mem_gb
        }
    }

    # Do kinship across all combinations of subsets
    Array[Int] split_ids = range(length(split_plink.bed_out))
    Array[Pair[Int, Int]] split_combos = cross(split_ids, split_ids)
    scatter(split_combo in split_combos){
        Int split_1 = split_combo.left
        Int split_2 = split_combo.right
        if(split_1 < split_2){
            call KING.kinship as pairwise_kinships{
                input:
                    bed_in = split_plink.bed_out[split_1],
                    bim_in = split_plink.bim_out[split_1],
                    fam_in = split_plink.fam_out[split_1],
                    bed_in_pair = split_plink.bed_out[split_2],
                    bim_in_pair = split_plink.bim_out[split_2],
                    fam_in_pair = split_plink.fam_out[split_2],
                    degree = degree,
                    output_basename = "${output_basename}.splitcombo.${split_1}.${split_2}",
                    cpu = king_split_cpu,
                    mem_gb = king_split_mem_gb
            }
        }
    }

    # Flatten to get all files in a single array
    # Have to do a select_all on pairwise combo outputs because we skipped identity combos
    call UTILS.flatten_string_array{
        input:
            array=[select_all(pairwise_kinships.kinship_output), subset_kinships.kinship_output]
    }

    output{
        Array[File] kinship_outputs = flatten_string_array.flat_array
    }

}