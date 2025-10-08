version 1.1

import "tsv_utils.wdl" as TSV
import "utils.wdl" as UTILS
import "generate_gwas_plots.wdl" as PLOT

workflow summarize_gwas_wf{

    input{
        File summary_stats_input
        String output_basename
        Float sample_maf_cutoff = 0.0
        Float pop_maf_cutoff = 0.0
        Float min_rsq = 0.0
        Float sig_alpha

        Int sample_maf_col = 7
        Int pop_maf_col = 8
        Int pvalue_col = 14
        Int rsq_col = 10

        String id_colname = "VARIANT_ID"
        String chr_colname = "CHR"
        String pos_colname = "POS"
        String p_colname = "P"

        String image_source = "docker"
        String? ecr_repo

        Int plot_mem_gb
    }

    # Format rsq
    call UTILS.format_float_sig_digits as format_rsq{
        input:
            input_float = min_rsq,
            significant_digits = 3,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Format sample maf
    call UTILS.format_float_sig_digits as format_sample_maf{
        input:
            input_float = sample_maf_cutoff,
            significant_digits = 3,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Format pop maf
    call UTILS.format_float_sig_digits as format_pop_maf{
        input:
            input_float = pop_maf_cutoff,
            significant_digits = 3,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Format pvalue
    call UTILS.format_float_sig_digits as format_pvalue{
        input:
            input_float = sig_alpha,
            significant_digits = 3,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Optionally fiter on Rsq if desired
    if (min_rsq > 0.0){
        call TSV.tsv_filter as filter_rsq{
            input:
                tsv_input = summary_stats_input,
                output_filename = output_basename + "_rsq_~{format_rsq.formatted_float}.tsv",
                filter_string = "--is-numeric '~{rsq_col}' --ge '~{rsq_col}:~{min_rsq}'",
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Optionally filter by sample MAF
    File rsq_summary_stats = select_first([filter_rsq.tsv_output, summary_stats_input])
    if(sample_maf_cutoff > 0.0){
        String sample_filter_string = "--is-numeric '~{sample_maf_col}' --ge '~{sample_maf_col}:~{sample_maf_cutoff}'"
        call TSV.tsv_filter as filter_sample_maf{
            input:
                tsv_input = rsq_summary_stats,
                output_filename = basename(rsq_summary_stats, ".tsv") + "_sample_maf_~{format_sample_maf.formatted_float}.tsv",
                filter_string = sample_filter_string,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Optionally filter by population MAF
    File sample_maf_summary_stats = select_first([filter_sample_maf.tsv_output, rsq_summary_stats])
    if(pop_maf_cutoff > 0.0){
        String pop_filter_string = "--is-numeric '~{pop_maf_col}' --ge '~{pop_maf_col}:~{pop_maf_cutoff}'"
        call TSV.tsv_filter as filter_pop_maf{
            input:
                tsv_input = sample_maf_summary_stats,
                output_filename = basename(sample_maf_summary_stats, ".tsv") + "_pop_MAF_~{format_pop_maf.formatted_float}.tsv",
                filter_string = pop_filter_string,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Check if any variants remain after filtering
    File summary_stats = select_first([filter_pop_maf.tsv_output, sample_maf_summary_stats])
    call UTILS.wc as count_variants{
        input:
            input_file = summary_stats,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int variant_count = if(count_variants.num_lines > 0) then count_variants.num_lines - 1 else 0

    # Plot Manhattan and QQ plots
    if (variant_count > 0){
        call PLOT.generate_gwas_plots as make_plots{
            input:
                summary_stats = summary_stats,
                col_id = id_colname,
                col_chromosome = chr_colname,
                col_position = pos_colname,
                col_p = p_colname,
                output_basename = basename(summary_stats, ".tsv"),
                mem_gb = plot_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Filter by pvalue to get a more interesting set of significant hits
    call TSV.tsv_filter as filter_pvalue{
        input:
            tsv_input = summary_stats,
            output_filename = basename(summary_stats, ".tsv") + "_p_~{format_pvalue.formatted_float}.tsv",
            filter_string = "--is-numeric '~{pvalue_col}' --le '~{pvalue_col}:~{sig_alpha}'",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Gzip summary stats for downstream easiness
    call UTILS.gzip as gzip_summary_stats{
        input:
            input_file = summary_stats,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.gzip as gzip_sig_summary_stats{
        input:
            input_file = filter_pvalue.tsv_output,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File summary_stats_output = gzip_summary_stats.output_file
        File sig_summary_stats_output = gzip_sig_summary_stats.output_file
        Array[File]? summary_plots = make_plots.plots
    }
}