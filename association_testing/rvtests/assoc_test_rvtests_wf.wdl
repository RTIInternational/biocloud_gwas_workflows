import "biocloud_gwas_workflows/association_testing/rvtests/rvtests_chr_wf" as RVCHR
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/generate_gwas_plots/generate_gwas_plots.wdl" as PLOT

workflow assoc_test_rvtests_wf{
    Array[File] vcfs_in
    Array[File] infos_in

    File pheno_file
    String pheno_name
    File covar_file
    Array[String] covars
    String dosage

    # For annotating with population MAF info
    Array[File] pop_maf_files
    String maf_population

    # Number of records per VCF split
    Int records_per_split = 50000

    # Set output basename
    String study_output_basename

    # Filtering cutoffs
    Float min_rsq

    # MAF Filter levels
    Array[Float]? user_maf_cutoffs
    Boolean filter_by_pop_maf
    Boolean filter_by_sample_maf

    # P-value cutoff for considering significant
    Float sig_alpha

    # RVTests options
    Boolean is_continuous
    Boolean inverseNormal = is_continuous
    Boolean useResidualAsPhenotype = is_continuous
    Boolean qtl = is_continuous

    Array[String] metaTestsMaybe = ["score"]
    Boolean? sex
    Boolean? multipleAllele
    String? xLabel

    # TODO: Make kinship matrices if necessary

    # Do RVTests chr workflow on each chromosome in parallel
    scatter(chr_index in range(length(vcfs_in))){
        call RVCHR.rvtests_chr_wf as rvtests{
            input:
                vcf_in = vcfs_in[chr_index],
                info_in = infos_in[chr_index],
                pheno_file = pheno_file,
                pheno_name = pheno_name,
                covar_file = covar_file,
                covars = covars,
                dosage = dosage,
                pop_maf_file = pop_maf_files[chr_index],
                maf_population = maf_population,
                records_per_split = records_per_split,
                metaTestsMaybe = metaTestsMaybe,
                sex = sex,
                multipleAllele = multipleAllele,
                xLabel = xLabel,
                inverseNormal = inverseNormal,
                useResidualAsPhenotype = useResidualAsPhenotype,
                qtl = qtl
        }
    }

    # Merge chr outputs into single file
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
            input_files = rvtests.summary_stats,
            output_dir_name = study_output_basename + ".rvtests_chr_output"
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = study_output_basename + ".rvtests.merged_assoc.tsv"
    }

    # Fiter on Rsq
    call TSV.tsv_filter as filter_rsq{
        input:
            tsv_input = cat_sumstats.tsv_output,
            output_filename = basename(cat_sumstats.tsv_output, ".tsv") + ".rsq_filter.tsv",
            filter_string = "--ge '10:${min_rsq}'"
    }


    # Plot QQ and Manhattan of GWAS hits filtered on imputation quality
    call PLOT.generate_gwas_plots as make_plots{
        input:
            summary_stats = filter_rsq.tsv_output,
            col_id = "VARIANT_ID",
            col_chromosome = "CHR",
            col_position = "POS",
            col_p = "P",
            output_basename = basename(filter_rsq.tsv_output, ".tsv")
        }

    # Filter by pvalue to get a more interesting set of significant hits
    call TSV.tsv_filter as filter_pvalue{
        input:
            tsv_input = filter_rsq.tsv_output,
            output_filename = basename(filter_rsq.tsv_output, ".tsv") + ".sighits.tsv",
            filter_string = "--lt '13:${sig_alpha}'"
    }

    # Make filtered datasets (Sample MAF, Pop MAF) and plot results for each
    if((defined(user_maf_cutoffs)) && (filter_by_pop_maf || filter_by_sample_maf)){
        scatter(maf_cutoff in user_maf_cutoffs){
            if((maf_cutoff > 0.0) && (maf_cutoff < 1.0)){
                String sample_filter = if (filter_by_sample_maf) then "--gt '6:${maf_cutoff}' " else ""
                String pop_filter = if (filter_by_pop_maf) then "--gt '8:${maf_cutoff}'" else ""
                String filter_string = sample_filter + pop_filter
                call TSV.tsv_filter as filter_maf{
                    input:
                        tsv_input = filter_rsq.tsv_output,
                        output_filename = basename(cat_sumstats.tsv_output, ".tsv") + "mafGT.${maf_cutoff}.tsv",
                        filter_string = filter_string
                }

                # Plot manhattan and QQ plots
                call PLOT.generate_gwas_plots as make_maf_plots{
                    input:
                        summary_stats = filter_maf.tsv_output,
                        col_id = "VARIANT_ID",
                        col_chromosome = "CHR",
                        col_position = "POS",
                        col_p = "P",
                        output_basename = basename(filter_maf.tsv_output, ".tsv")
                }

                # Filter by pvalue to get a more interesting set of significant hits
                call TSV.tsv_filter as filter_maf_pvalue{
                    input:
                        tsv_input = filter_maf.tsv_output,
                        output_filename = basename(filter_maf.tsv_output, ".tsv") + ".sighits.tsv",
                        filter_string = "--lt '13:${sig_alpha}'"
                }
            }
        }
    }

    output{
        File rsq_filtered_summary_stats = filter_rsq.tsv_output
        File sig_rsq_filtered_summary_stats = filter_pvalue.tsv_output
        Array[File] rsq_filtered_gwas_plots = make_plots.plots
        Array[File] maf_summary_stats = select_all(filter_maf.tsv_output)
        Array[File] sig_maf_summary_stats = select_all(filter_maf_pvalue.tsv_output)
        Array[Array[File]] maf_gwas_plots = select_all(make_maf_plots.plots)

    }
}