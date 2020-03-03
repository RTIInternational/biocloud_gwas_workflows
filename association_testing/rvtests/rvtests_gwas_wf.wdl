import "biocloud_gwas_workflows/association_testing/rvtests/rvtests_gwas_chr_wf.wdl" as RVCHR
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/summarize_gwas_wf.wdl" as SUM

workflow rvtests_gwas_wf{
    Array[File] vcfs_in
    Array[File] infos_in

    File pheno_file
    String pheno_name
    File? covar_file
    Array[String]? covars
    String dosage

    # For annotating with population MAF info
    Array[File] pop_maf_files
    String maf_population

    # Number of records per VCF split
    Int records_per_split = 100000

    # Set output basename
    String study_output_basename

    # Filtering cutoffs
    Float min_rsq

    # MAF Filter levels
    Array[Float]? maf_cutoffs
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
        call RVCHR.rvtests_gwas_chr_wf as rvtests{
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
            output_dir_name = study_output_basename + "_rvtests_chr_output"
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = study_output_basename + ".rvtests.MetaAssoc.tsv"
    }

    # Fiter on Rsq if desired
    if (min_rsq > 0.0){
        call TSV.tsv_filter as filter_rsq{
            input:
                tsv_input = cat_sumstats.tsv_output,
                output_filename = basename(cat_sumstats.tsv_output, ".tsv") + ".rsq.tsv",
                filter_string = "--ge '10:${min_rsq}'"
        }
    }

    # Make manhattan, QQ plots of sumstats with no MAF filtering
    File sumstats_file = select_first([filter_rsq.tsv_output, cat_sumstats.tsv_output])
    call SUM.summarize_gwas_wf as summarize_rsq_sumstats{
        input:
            summary_stats_input = sumstats_file,
            output_basename = basename(sumstats_file, ".tsv"),
            sig_alpha = sig_alpha
    }

    # Summarize results after applying user-defined MAF cutoffs
    if((defined(maf_cutoffs)) && (filter_by_pop_maf || filter_by_sample_maf)){
        Array[Float] actual_cutoffs = select_first([maf_cutoffs, []])
        scatter(maf_cutoff in actual_cutoffs){
            Float pop_maf_cutoff = if filter_by_pop_maf then maf_cutoff else 0.0
            Float sample_maf_cutoff = if filter_by_sample_maf then maf_cutoff else 0.0
            call SUM.summarize_gwas_wf as summarize_maf_sumstats{
                input:
                    summary_stats_input = sumstats_file,
                    output_basename = basename(sumstats_file, ".tsv"),
                    sig_alpha = sig_alpha,
                    sample_maf_cutoff = sample_maf_cutoff,
                    pop_maf_cutoff = pop_maf_cutoff
            }
        }
    }

    output{
        File summary_stats = sumstats_file
        File sig_summary_stats = summarize_rsq_sumstats.sig_summary_stats_output
        Array[File] summary_stats_plots = summarize_rsq_sumstats.summary_plots
        Array[File]? maf_summary_stats = summarize_maf_sumstats.summary_stats_output
        Array[File]? sig_maf_summary_stats = summarize_maf_sumstats.sig_summary_stats_output
        Array[Array[File]]? maf_summary_stats_plots = summarize_maf_sumstats.summary_plots
    }
}