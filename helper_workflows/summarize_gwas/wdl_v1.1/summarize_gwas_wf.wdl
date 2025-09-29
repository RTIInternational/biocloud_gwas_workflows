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

    # Optionally fiter on Rsq if desired
    if (min_rsq > 0.0){
        call TSV.tsv_filter as filter_rsq{
            input:
                tsv_input = summary_stats_input,
                output_filename = output_basename + ".rsq.~{min_rsq}.tsv",
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
                output_filename = basename(rsq_summary_stats, ".tsv") + ".sampleMAF.~{sample_maf_cutoff}.tsv",
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
                output_filename = basename(sample_maf_summary_stats, ".tsv") + ".popMAF.~{pop_maf_cutoff}.tsv",
                filter_string = pop_filter_string,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Plot Manhattan and QQ plots
    File summary_stats = select_first([filter_pop_maf.tsv_output, sample_maf_summary_stats])
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

    # Filter by pvalue to get a more interesting set of significant hits
    call TSV.tsv_filter as filter_pvalue{
        input:
            tsv_input = summary_stats,
            output_filename = basename(summary_stats, ".tsv") + ".p.~{sig_alpha}.tsv",
            filter_string = "--is-numeric '~{pvalue_col}' --le '~{pvalue_col}:~{sig_alpha}'",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Gzip summary stats for downstream easiness
    call UTILS.gzip{
        input:
            input_file = summary_stats,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File summary_stats_output = summary_stats
        File gzipped_summary_stats_output = gzip.output_file
        File sig_summary_stats_output = filter_pvalue.tsv_output
        Array[File] summary_plots = make_plots.plots
    }
}