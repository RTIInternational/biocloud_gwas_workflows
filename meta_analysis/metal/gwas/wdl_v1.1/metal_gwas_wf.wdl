version 1.1

import "metal.wdl" as METAL
import "rti-tsv-utils.wdl" as JOIN
import "tsv_utils.wdl" as TSV
import "generate_gwas_plots.wdl" as PLOT

workflow metal_gwas_wf {

    input{
        # METAL input file parameters
        Array[File] sum_stats_files
        Array[String] separators
        Array[String] marker_col_names
        Array[String] chrom_col_names
        Array[String] pos_col_names
        Array[String] ref_allele_col_names
        Array[String] alt_allele_col_names
        Array[String] effect_col_names
        Array[String] freq_col_names = []

        # METAL input file sample size weighted meta parameters
        Array[String] pvalue_col_names = []
        Array[String] weight_col_names = []

        # METAL input file inverse variance weighted meta parameters
        Array[String] std_err_col_names = []

        # METAL analysis parameters
        String metal_out_prefix
        String metal_out_suffix = "metal"
        String scheme = "STDERR"
        String ?analyze     # Set to "HETEROGENEITY" to run heterogeneity analysis; otherwise, leave unset for standard meta-analysis
        String ?genomic_control     # Set to "ON" to apply genomic control; otherwise, leave unset

        # Other METAL parameters
        String column_counting = "STRICT"
        String average_freq = "ON"
        String min_max_freq = "ON"
        String track_positions = "ON"

        # Parameters for joining METAL results with sum stats files
        Array[String] sum_stats_files_labels
        Int metal_sum_stats_join_chunk_size = 2000000
        Int metal_sum_stats_join_cpu = 1
        Int metal_sum_stats_join_mem_gb = 2

        # Parameters for p-value filter
        Float p_threshold = 0.001
        Int p_filter_cpu = 1
        Int p_filter_mem_gb = 2

        # METAL resources
        Int metal_cpu = 1
        Int metal_mem_gb = 2

        # Plot resources
        Int plot_cpu = 1
        Int plot_mem_gb = 2

        # Runtime options
        String image_source = "docker"
        String? ecr_repo
    }

    # Create METAL command file
    call METAL.metal as metal {
        input:
            sum_stats_files = sum_stats_files,
            separators = separators,
            marker_col_names = marker_col_names,
            chrom_col_names = chrom_col_names,
            pos_col_names = pos_col_names,
            ref_allele_col_names = ref_allele_col_names,
            alt_allele_col_names = alt_allele_col_names,
            effect_col_names = effect_col_names,
            freq_col_names = freq_col_names,
            pvalue_col_names = pvalue_col_names,
            weight_col_names = weight_col_names,
            std_err_col_names = std_err_col_names,
            metal_out_prefix = metal_out_prefix,
            metal_out_suffix = metal_out_suffix,
            scheme = scheme,
            genomic_control = genomic_control,
            analyze = analyze,
            column_counting = column_counting,
            average_freq = average_freq,
            min_max_freq = min_max_freq,
            track_positions = track_positions,
            cpu = metal_cpu,
            mem_gb = metal_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Create input file for plots
    call TSV.tsv_select as plot_tsv {
        input:
            tsv_input = metal.metal_results,
            output_filename = "plot.tsv",
            fields = ["MarkerName", "Chromosome", "Position", "P"],
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Generate plots
    call PLOT.generate_gwas_plots as generate_plots{
        input:
            summary_stats = plot_tsv.tsv_output,
            col_id = "MarkerName",
            col_chromosome = "Chromosome",
            col_position = "Position",
            col_p = "P",
            output_basename = "~{metal_out_prefix}",
            cpu = plot_cpu,
            mem_gb = plot_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get right remove cols for sum stats files
    call combine_array_elements as sum_stats_remove_cols {
        input:
            array1 = chrom_col_names,
            array2 = pos_col_names,
            separator = ",",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Append dataset-specific sum stats to METAL results
    call JOIN.tsv_join as metal_sum_stats_join {
        input:
            left_file = metal.metal_results,
            left_on = "MarkerName",
            right_files = sum_stats_files,
            right_ons = marker_col_names,
            right_suffixes = sum_stats_files_labels,
            right_seps = separators,
            right_remove_cols = sum_stats_remove_cols.output_array,
            hows = ["left"],
            out_prefix = metal_out_prefix + "_with_all_sum_stats",
            sort = false,
            chunk_size = metal_sum_stats_join_chunk_size,
            mem_gb = metal_sum_stats_join_mem_gb,
            cpu = metal_sum_stats_join_cpu,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Filter by p-value
    call TSV.tsv_filter as p_value_filter {
        input:
            tsv_input = metal_sum_stats_join.out_tsv,
            output_filename = metal_out_prefix + "_p_lt_" + p_threshold + ".tsv",
            filter_string = "--lt P:" + p_threshold,
            cpu = p_filter_cpu,
            mem_gb = p_filter_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File meta_sum_stats = metal.metal_results
        File meta_with_full_sum_stats = metal_sum_stats_join.out_tsv
        File meta_with_full_sum_stats_p_filter = p_value_filter.tsv_output
        File metal_info = metal.metal_info
        File metal_log = metal.metal_log
        Array[File] plots = generate_plots.plots
    }

}

task combine_array_elements {
    input {
        Array[String] array1
        Array[String] array2
        String separator = ","
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 2
    }

    command <<<
        set -e
        array1_str="~{sep(' ', array1)}"
        array1=($array1_str)
        array1_length=${#array1[@]}
        array2_str="~{sep(' ', array2)}"
        array2=($array2_str)
        array2_length=${#array2[@]}
        if [ $array1_length -ne $array2_length ]; then
            echo "Error: Input arrays must have the same length" >&2
            exit 1
        fi
        for (( i=0; i<$array1_length; i++)); do
            echo "${array1[$i]}~{separator}${array2[$i]}"
        done > combined.tsv
    >>>

    output {
        Array[String] output_array = read_lines("combined.tsv")
    }
    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }
}