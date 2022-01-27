import "biocloud_wdl_tools/make_gwas_summary_stats/make_gwas_summary_stats.wdl" as SUMSTAT

workflow test_make_gwas_summary_stats{
    File rvtests_sumstats_file
    File gem_sumstats_file
    File pop_maf_file
    File info_file
    File mfi_file

    call SUMSTAT.make_gwas_summary_stats as rvtests_anno_sumstats_eur{
        input:
            file_in_summary_stats = rvtests_sumstats_file,
            file_in_info = info_file,
            file_in_pop_mafs = pop_maf_file,
            file_in_summary_stats_format = "rvtests",
            population = 'eur',
            file_out_prefix = "rvtests_chr22_sumstats_eur"
    }

    call SUMSTAT.make_gwas_summary_stats as gem_anno_sumstats_eur{
        input:
            file_in_summary_stats = gem_sumstats_file,
            file_in_info = mfi_file,
            file_in_summary_stats_format = "gem",
            file_in_info_format = "mfi",
            file_out_prefix = "gem_chr22_sumstats_eur"
    }

    output{
        File rvtests_annotated_sumstats_eur = rvtests_anno_sumstats_eur.output_file
        File gem_annotated_sumstats_eur = gem_anno_sumstats_eur.output_file
    }
}