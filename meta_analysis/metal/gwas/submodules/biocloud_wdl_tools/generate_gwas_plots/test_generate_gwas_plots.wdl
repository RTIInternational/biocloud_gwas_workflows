import "biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT

workflow test_generate_gwas_plots{
    File rvtests_sumstats_file
    File probabel_sumstats_file
    Int max_retries = 1

    call PLOT.generate_gwas_plots as plot_rvtests{
        input:
            summary_stats = rvtests_sumstats_file,
            col_id = "FAKE_MARKER_NAME",
            col_chromosome = "CHROM",
            col_position = "POS",
            col_p = "PVALUE",
            output_basename = "rvtests",
            max_retries = max_retries

    }

    call PLOT.generate_gwas_plots as plot_probabel{
        input:
            summary_stats = probabel_sumstats_file,
            col_id = "SNP",
            col_chromosome = "chr",
            col_position = "pos",
            col_p = "pval",
            output_basename = "probabel",
            max_retries = max_retries
    }
    output{
        Array[File] rvtests_plots = plot_rvtests.plots
        Array[File] probabel_plots = plot_probabel.plots
    }
}

