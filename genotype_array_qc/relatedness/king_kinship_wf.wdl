import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/king/king.wdl" as KING
import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

workflow king_kinship_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    # Degree of relatedness cutoff for kinship file
    Int? degree

    # Number of subsets to make for splitting
    Int num_splits = 4

    Int plink_cpu = 1
    Int plink_mem_gb = 2

    Int king_split_cpu = 2
    Int king_split_mem_gb = 4
    
    String docker_ubuntu = "404545384114.dkr.ecr.us-east-1.amazonaws.com/ubuntu:18.04"
    String docker_pigz = "404545384114.dkr.ecr.us-east-1.amazonaws.com/rtibiocloud/pigz:v2.4_b243f9"
    String docker_plink2_0 = "404545384114.dkr.ecr.us-east-1.amazonaws.com/rtibiocloud/plink:v2.0_4d3bad3"
    String docker_king = "404545384114.dkr.ecr.us-east-1.amazonaws.com/rtibiocloud/king:v2.24_9b4c1b9"

    # Split plink sample file into however many chunks
    call SPLIT.split_file as split_fam{
        input:
            input_file = fam_in,
            output_basename = basename(fam_in, ".fam"),
            output_extension = ".fam",
            num_splits = num_splits,
            docker = docker_pigz
    }

    # Create split plink files
    scatter(split_index in range(length(split_fam.output_files))){

        # Generate keep file from split fam
        call UTILS.cut{
            input:
                input_file = split_fam.output_files[split_index],
                output_filename = "${output_basename}.split.${split_index}.keep",
                args = "-f 1,2",
                docker = docker_ubuntu
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
                mem_gb = plink_mem_gb,
                docker = docker_plink2_0
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
                mem_gb = king_split_mem_gb,
                docker = docker_king
        }
    }

    # Compute kinship between all pairwise combinations of sample subsets if parallelizing
    if(num_splits > 1){
        # Make cross product of all possible pairwise combinations
        Array[Int] split_ids = range(length(split_plink.bed_out))
        Array[Pair[Int, Int]] split_combos = cross(split_ids, split_ids)
        scatter(split_combo in split_combos){
            Int split_1 = split_combo.left
            Int split_2 = split_combo.right
            # Neat trick to only do unique combos
            # E.g. if you cross [1,2] with [1,2] split_combos will contain (1,2) and (2,1)
            # This automatically eliminates identity and duplicate pairwise runs
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
                        mem_gb = king_split_mem_gb,
                        docker = docker_king
                }
            }
        }

        # Flatten to get all files in a 1-D array
        call UTILS.flatten_string_array{
            input:
                array=[select_all(pairwise_kinships.kinship_output), subset_kinships.kinship_output],
                docker = docker_ubuntu
        }

        # Zip into single tarball for tsv-concat
        call COLLECT.collect_large_file_list_wf as collect_kinships{
            input:
                input_files = flatten_string_array.flat_array,
                output_dir_name = "${output_basename}_kinships"
        }

        # Concat all kinship files (preserving header) into single kinship file
        call TSV.tsv_append as cat_kinships{
            input:
                tsv_inputs_tarball = collect_kinships.output_dir,
                output_filename = "${output_basename}.merged.kinship.kin0"
        }
    }

    File king_kinship = select_first([cat_kinships.tsv_output, subset_kinships.kinship_output[0]])
    # Prune resulting kinship file to generate minimal set of related individuals
    call KING.prune_related_samples{
        input:
            kinship_in = king_kinship,
            output_basename = "${output_basename}.pruned",
            mem_gb = ceil(size(king_kinship, "GB")) + 1
    }

    output{
        # Select either the merged kinship file (num_splits > 1) or the first shard (num_splits == 1)
        File raw_kinship_output = king_kinship
        File related_samples = prune_related_samples.related_samples
        File annotated_kinship_output = prune_related_samples.annotated_kinship_out
    }

}
