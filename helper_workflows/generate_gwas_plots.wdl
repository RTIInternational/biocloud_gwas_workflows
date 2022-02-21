import "biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT

workflow generate_gwas_plots{

    File sumstats_file
    String col_id
    String col_chromosome
    String col_position
    String col_p
    String? col_variant_type
    String output_basename

    Boolean? in_header = true
    Boolean? in_csv
    File? highlight_list

    Boolean? generate_manhattan_plot = true
    Boolean? generate_qq_plot = true
    Boolean? generate_snp_manhattan_plot
    Boolean? generate_indel_manhattan_plot
    Boolean? generate_snp_indel_manhattan_plot
    Boolean? generate_snp_qq_plot
    Boolean? generate_indel_qq_plot
    Boolean? generate_snp_indel_qq_plot

    Boolean? qq_lambda = true
    Boolean? qq_lines = true
    Boolean? qq_significance_line
    Int? manhattan_ylim
    Boolean? manhattan_no_line
    String? manhattan_odd_chr_color
    String? manhattan_even_chr_color
    String? manhattan_highlight_chr_color
    Float? manhattan_significance_value

    # Runtime options
    Int? cpu
    Int? mem_gb


    Int max_retries = 3

    call PLOT.generate_gwas_plots as plot{
        input:
            summary_stats = sumstats_file,
            col_id = col_id,
            col_chromosome = col_chromosome,
            col_position = col_position,
            col_p = col_p,
            variant_type_colname = col_variant_type,
            output_basename = output_basename,
            in_header = in_header,
            in_csv = in_csv,
            highlight_list = highlight_list,
            generate_manhattan_plot = generate_manhattan_plot,
            generate_qq_plot = generate_qq_plot,
            generate_snp_manhattan_plot = generate_snp_manhattan_plot,
            generate_indel_manhattan_plot = generate_indel_manhattan_plot,
            generate_snp_indel_manhattan_plot = generate_snp_indel_manhattan_plot,
            generate_snp_qq_plot = generate_snp_qq_plot,
            generate_indel_qq_plot = generate_indel_qq_plot,
            generate_snp_indel_qq_plot = generate_snp_indel_qq_plot,
            qq_lambda = qq_lambda,
            qq_lines = qq_lines,
            qq_significance_line = qq_significance_line,
            manhattan_ylim = manhattan_ylim,
            manhattan_no_line = manhattan_no_line,
            manhattan_odd_chr_color = manhattan_odd_chr_color,
            manhattan_even_chr_color = manhattan_even_chr_color,
            manhattan_highlight_chr_color = manhattan_highlight_chr_color,
            manhattan_significance_value = manhattan_significance_value,
            cpu = cpu,
            mem_gb = mem_gb
    }

    output{
        Array[File] plots = plot.plots
    }
}

