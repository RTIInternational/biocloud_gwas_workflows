import "biocloud_gwas_workflows/association_testing/gem/gem_gwas_chr_wf.wdl" as GEM_CHR
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/helper_workflows/summarize_gwas_wf_v2.wdl" as SUM

workflow gem_gwas_bgen_wf{

    Array[String] chrs

    # GEM Input/Output File Options:
    File pheno_file
    Array[File] in_bgen_files
    Array[File] in_sample_files
    String out_file_prefix
    String output_style = "minimum"

    # GEM Phenotype File Options:
    String sampleid_name
    String pheno_name
    String? exposure_names
    String? int_covar_names
    String? covar_names
    Int robust = 0
    Float? tol
    String? delim
    String? missing_value
    Int? center
    Int? scale
    String? categorical_names
    Int? cat_threshold

    # GEM Filtering Options:
    Float? maf
    Float? miss_geno_cutoff
    String? include_snp_file

    # GEM Performance Options:
    Int? gem_threads
    Int? stream_snps

    # GEM Runtime attributes
    Int gem_cpu = 1
    Int gem_mem_gb = 2

    # Options for make_gwas_summary_stats
    Array[File] info_files
    String info_file_format

    # ID conversion options
    Array[File] id_ref_files
    String in_missing_allele
    String in_deletion_allele
    String ref_deletion_allele = "."
    String id_label = ""
    Int id_conversion_cpu = 1
    Int id_conversion_mem_gb = 1

    # MAF filtering options
    Boolean filter_by_sample_maf        # Boolean switch for whether to filter by sample MAF
    Boolean filter_by_pop_maf           # Boolean switch for whether to filter by population MAF
    Array[Float]? maf_cutoffs           # Optional MAF cutoffs levels used to generate summary stats/graphs across a range of MAF filters

    # Imputation quality filtering options
    Boolean filter_by_rsq           # Boolean switch for whether to merge RSQ into results and whether to filter by rsq
    Array[Float]? rsq_cutoffs       # Imputation quality filtering cutoffs. Final outputs will have R2 >= this value

    # P-value filtering
    Boolean filter_by_pvalue            # Boolean switch for whether to filter by p-value
    Float sig_alpha                     # P-value threshold for filtering
    Array[String] non_robust_p_value_fields = if(output_style == "minimum" && robust == 1) then [] else ["P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint"]
    Array[String] robust_p_value_fields = if(robust == 1) then ["robust_P_Value_Marginal", "robust_P_Value_Interaction", "robust_P_Value_Joint"] else []
    Array[String] p_value_fields = flatten([non_robust_p_value_fields, robust_p_value_fields])

    # Mem for plotting
    Int plot_mem_gb

    # Do GEM chr workflow on each chromosome in parallel
    scatter(index in range(length(in_bgen_files))){

        call GEM_CHR.gem_gwas_chr_wf as GEM{
            input:
                pheno_file = pheno_file,
                out = out_file_prefix + '_chr' + chrs[index],
                bgen = in_bgen_files[index],
                sample = in_sample_files[index],
                output_style = output_style,
                sampleid_name = sampleid_name,
                pheno_name = pheno_name,
                exposure_names = exposure_names,
                int_covar_names = int_covar_names,
                covar_names = covar_names,
                robust = robust,
                tol = tol,
                delim = delim,
                missing_value = missing_value,
                center = center,
                scale = scale,
                categorical_names = categorical_names,
                cat_threshold = cat_threshold,
                maf = maf,
                miss_geno_cutoff = miss_geno_cutoff,
                include_snp_file = include_snp_file,
                threads = gem_threads,
                stream_snps = stream_snps,
                gem_cpu = gem_cpu,
                gem_mem_gb = gem_mem_gb,
                info_file = info_files[index],
                info_file_format = info_file_format
        }

        # Convert ids in summary stats output to standard ids
        call IDCONVERT.convert_variant_ids{
            input:
                chr = chrs[index],
                in_file = GEM.summary_stats,
                in_header = 1,
                in_sep = "tab",
                ref = id_ref_files[index],
                in_id_col = 0,
                in_chr_col = 1,
                in_pos_col = 2,
                in_a1_col = 3,
                in_a2_col = 4,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = ref_deletion_allele,
                output_filename = out_file_prefix + '_chr' + chrs[index] + '_' + id_label + ".tsv.gz",
                output_compression = "gzip",
                cpu = id_conversion_cpu,
                mem_gb = id_conversion_mem_gb
        }

    }

    # Merge chr outputs into single file
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
            input_files = convert_variant_ids.output_file,
            output_dir_name = out_file_prefix + "_gem_chr_output"
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = out_file_prefix + ".tsv"
    }

    # Generate summaries for each possible combo of MAF/Rsq cutoff
    Array[Float] actual_maf_cutoffs = select_first([maf_cutoffs, [0.0]])
    Array[Float] actual_rsq_cutoffs = select_first([rsq_cutoffs, [0.0]])
    Array[Pair[Float, Float]] filter_combos = cross(actual_maf_cutoffs, actual_rsq_cutoffs)
    scatter(filter_combo in filter_combos){

        # Get the MAF/Rsq cutoffs for this combo
        Float min_rsq = filter_combo.right
        Float min_maf = filter_combo.left

        # Set the pop/sample mafs based on whether we want to filter on sample or population MAF
        Float pop_maf_cutoff = if filter_by_pop_maf then min_maf else 0.0
        Float sample_maf_cutoff = if filter_by_sample_maf then min_maf else 0.0

        # Generate summary of
        call SUM.summarize_gwas_wf as summarize_filtered_sumstats{
            input:
                summary_stats_input = cat_sumstats.tsv_output,
                output_basename = basename(cat_sumstats.tsv_output, ".tsv"),
                min_rsq = min_rsq,
                sample_maf_cutoff = sample_maf_cutoff,
                pop_maf_cutoff = pop_maf_cutoff,
                sig_alpha = sig_alpha,
                pvalue_colnames = p_value_fields,
                plot_mem_gb = plot_mem_gb
        }

    }

    output{
        File raw_summary_stats = cat_sumstats.tsv_output
        Array[File] filtered_summary_stats = summarize_filtered_sumstats.summary_stats_output
        Array[File] gzipped_filtered_summary_stats = summarize_filtered_sumstats.gzipped_summary_stats_output
        Array[Array[File]] sig_filtered_summary_stats = summarize_filtered_sumstats.sig_summary_stats_output
        Array[Array[Array[File]]] filtered_summary_stats_plots = summarize_filtered_sumstats.summary_plots
    }
}