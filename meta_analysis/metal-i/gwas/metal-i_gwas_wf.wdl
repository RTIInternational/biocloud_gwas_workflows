import "biocloud_gwas_workflows/biocloud_wdl_tools/metal-i/metal-i.wdl" as METAL
import "biocloud_gwas_workflows/biocloud_wdl_tools/rti-tsv-utils/rti-tsv-utils.wdl" as JOIN
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT

workflow metal_i_gwas_wf {

    # METAL input file parameters
    Array[File] sum_stats_files
    String separators = "WHITESPACE"
    String marker_col_names = "VARIANT_ID"
    String ref_allele_col_names = "REF"
    String alt_allele_col_names = "ALT"
    String effect_col_names = "ALT_EFFECT"
    String freq_col_names = "ALT_AF"

    # METAL input file sample size weighted meta parameters
    String pvalue_col_names = "P"
    String weight_col_names = "N"

    # METAL input file inverse variance weighted meta parameters
    String std_err_col_names = "SE"

    # METAL input file interaction parameters
    String ?int_effect_col_names
    String ?int_std_err_col_names
    String ?int_cov_col_names

    # METAL analysis parameters
    String metal_out_prefix
    String ?metal_out_suffix = "metal"
    String scheme = "INTERACTION"
    String ?analyze
    String ?genomic_control

    # Other METAL parameters
    String ?column_counting
    String average_freq = "ON"
    String min_max_freq = "ON"

    # METAL resources
    Int metal_cpu = 1
    Int metal_mem_gb = 2

    # Parameters for joining METAL results with chr & pos
    String variant_chr_pos_file
    String variant_chr_pos_file_id_col = "VARIANT_ID"
    String variant_chr_pos_file_chr_col = "CHR"
    String variant_chr_pos_file_pos_col = "POS"
    String variant_chr_pos_file_sep = "whitespace"
    Int metal_chr_pos_join_chunk_size = 5000000
    Int metal_chr_pos_join_mem_gb = 4

    # Parameters for joining METAL results with sum stats files
    String sum_stats_files_labels
    Int metal_sum_stats_join_chunk_size = 2000000
    Int metal_sum_stats_join_cpu = 1
    Int metal_sum_stats_join_mem_gb = 2

    # Plot resources
    Int plot_cpu = 1
    Int plot_mem_gb = 2

    # Parameters for p-value filter
    Float p_threshold = 0.001
    Int p_filter_cpu = 1
    Int p_filter_mem_gb = 2

    # Create METAL command file
    call METAL.metal_i as metal_i {
        input:
            sum_stats_files = sum_stats_files,
            separators = separators,
            marker_col_names = marker_col_names,
            ref_allele_col_names = ref_allele_col_names,
            alt_allele_col_names = alt_allele_col_names,
            effect_col_names = effect_col_names,
            freq_col_names = freq_col_names,
            pvalue_col_names = pvalue_col_names,
            weight_col_names = weight_col_names,
            std_err_col_names = std_err_col_names,
            int_effect_col_names = int_effect_col_names,
            int_std_err_col_names = int_std_err_col_names,
            int_cov_col_names = int_cov_col_names,
            metal_out_prefix = metal_out_prefix,
            metal_out_suffix = metal_out_suffix,
            scheme = scheme,
            genomic_control = genomic_control,
            analyze = analyze,
            column_counting = column_counting,
            average_freq = average_freq,
            min_max_freq = min_max_freq,
            cpu = metal_cpu,
            mem_gb = metal_mem_gb
    }

    # Join METAL results with chr & pos
    call JOIN.tsv_join as metal_chr_pos_join {
        input:
            left_file = variant_chr_pos_file,
            left_on = variant_chr_pos_file_id_col,
            left_sep = variant_chr_pos_file_sep,
            left_cols = variant_chr_pos_file_id_col + "," + variant_chr_pos_file_chr_col + "," + variant_chr_pos_file_pos_col,
            right_files = [metal_i.metal_results],
            right_ons = "MarkerName",
            right_seps = "whitespace",
            right_suffixes = "_meta",
            hows = "inner",
            out_prefix = metal_out_prefix,
            sort = false,
            chunk_size = 5000000,
            mem_gb = metal_chr_pos_join_mem_gb
    }

    # Create input file for plots
    call TSV.tsv_select as plot_tsv {
        input:
            tsv_input = metal_chr_pos_join.out_tsv,
            output_filename = "plot.tsv",
            fields = [variant_chr_pos_file_id_col, variant_chr_pos_file_chr_col, variant_chr_pos_file_pos_col, "P_meta"]
    }

    # Generate plots
    call PLOT.generate_gwas_plots as generate_plots{
        input:
            summary_stats = plot_tsv.tsv_output,
            col_id = variant_chr_pos_file_id_col,
            col_chromosome = variant_chr_pos_file_chr_col,
            col_position = variant_chr_pos_file_pos_col,
            col_p = "P_meta",
            output_basename = "${metal_out_prefix}",
            mem_gb = plot_mem_gb
    }

    # Append dataset-specific sum stats to METAL results
    call JOIN.tsv_join as metal_sum_stats_join {
        input:
            left_file = metal_chr_pos_join.out_tsv,
            left_on = variant_chr_pos_file_id_col,
            right_files = sum_stats_files,
            right_ons = marker_col_names,
            right_seps = separators,
            right_suffixes = sum_stats_files_labels,
            hows = "left",
            out_prefix = metal_out_prefix + "_with_all_sum_stats",
            sort = false,
            chunk_size = 2000000,
            mem_gb = metal_sum_stats_join_mem_gb
    }

    # Filter by p-value
    call TSV.tsv_filter as p_value_filter {
        input:
            tsv_input = metal_sum_stats_join.out_tsv,
            output_filename = metal_out_prefix + "_p_lt_" + p_threshold + ".tsv",
            filter_string = "--lt P_meta:" + p_threshold,
            cpu = p_filter_cpu,
            mem_gb = p_filter_mem_gb
    }

    output{
        File meta_sum_stats = metal_chr_pos_join.out_tsv
        File meta_with_full_sum_stats = metal_sum_stats_join.out_tsv
        File meta_with_full_sum_stats_p_filter = p_value_filter.tsv_output
        File metal_info = metal_i.metal_info
        File metal_log = metal_i.metal_log
        Array[File] plots = generate_plots.plots
    }

}