import "biocloud_gwas_workflows/association_testing/genesis/genesis_gwas_chr_wf.wdl" as GENESIS_CHR
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/helper_workflows/summarize_gwas_wf_v2.wdl" as SUM

workflow genesis_gwas_wf{

    String container_source = "docker"

    Array[File] genotype_files
    Array[File] variant_lists = []
    Boolean use_variant_lists = false
    File output_file_prefix
    Array[String] chrs
    File pheno_file
    String pheno_name
    # File? covar_file
    Array[String]? covars
    String geno_format
    String family
    String? gxe

    Int genesis_cpu = 1
    Int genesis_mem_gb = 1

    # Split options
    Boolean split_gds = false
    Int? chunk_size
    String? variant_id_field
    Int? split_by_variant_cpu
    Int? split_by_variant_mem_gb

    # ID conversion parameters
    Array[File] id_ref_files
    String in_missing_allele
    String in_deletion_allele
    String ref_deletion_allele = "."
    String id_label = ""
    Int id_conversion_cpu = 1
    Int id_conversion_mem_gb = 1

    # TSV utils options
    Int tsv_append_mem_gb = 4
    
    # Options for make_gwas_summary_stats
    Array[File] info_files
    String info_file_format

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
    Array[String] p_value_fields = ["Score.pval"]

    # Mem for plotting
    Int plot_mem_gb

    # Do genesis chr workflow on each chromosome in parallel
    scatter(index in range(length(genotype_files))){

        if (!use_variant_lists) {
            call GENESIS_CHR.genesis_gwas_chr_wf as genesis{
                input:
                    container_source = container_source,
                    file_in_geno = genotype_files[index],
                    file_out_prefix = output_file_prefix + '_chr' + chrs[index],
                    file_in_pheno = pheno_file,
                    geno_format = geno_format,
                    pheno = pheno_name,
                    covars = covars,
                    family = family,
                    gxe = gxe,
                    chr = chrs[index],
                    genesis_cpu = genesis_cpu,
                    genesis_mem_gb = genesis_mem_gb,
                    split_gds = split_gds,
                    chunk_size = chunk_size,
                    variant_id_field = variant_id_field,
                    file_in_info = info_files[index],
                    file_in_info_format = info_file_format,
                    split_by_variant_cpu = split_by_variant_cpu,
                    split_by_variant_mem_gb = split_by_variant_mem_gb,
                    tsv_append_mem_gb = tsv_append_mem_gb
            }
        }

        if (use_variant_lists) {
            call GENESIS_CHR.genesis_gwas_chr_wf as genesis_variant_list{
                input:
                    container_source = container_source,
                    file_in_geno = genotype_files[index],
                    file_in_variant_list = variant_lists[index],
                    use_variant_list = use_variant_lists,
                    file_out_prefix = output_file_prefix + '_chr' + chrs[index],
                    file_in_pheno = pheno_file,
                    geno_format = geno_format,
                    pheno = pheno_name,
                    covars = covars,
                    family = family,
                    gxe = gxe,
                    chr = chrs[index],
                    genesis_cpu = genesis_cpu,
                    genesis_mem_gb = genesis_mem_gb,
                    split_gds = split_gds,
                    chunk_size = chunk_size,
                    variant_id_field = variant_id_field,
                    file_in_info = info_files[index],
                    file_in_info_format = info_file_format,
                    split_by_variant_cpu = split_by_variant_cpu,
                    split_by_variant_mem_gb = split_by_variant_mem_gb,
                    tsv_append_mem_gb = tsv_append_mem_gb
            }
        }

        # Convert ids in summary stats output to standard ids
        File summary_stats = select_first([genesis.summary_stats, genesis_variant_list.summary_stats])
        call IDCONVERT.convert_variant_ids{
            input:
                chr = chrs[index],
                in_file = summary_stats,
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
                output_filename = basename(summary_stats, ".tsv") + id_label + ".tsv.gz",
                output_compression = "gzip",
                container_source = container_source,
                cpu = id_conversion_cpu,
                mem_gb = id_conversion_mem_gb
        }

    }

    # Merge chr outputs into single file
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
            input_files = convert_variant_ids.output_file,
            output_dir_name = output_file_prefix + "_genesis_chr_output",
            container_source = container_source
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = output_file_prefix + ".tsv",
            container_source = container_source
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
                plot_mem_gb = plot_mem_gb,
                container_source = container_source
        }

    }

    output{
        File unfiltered_sumstats = cat_sumstats.tsv_output
        Array[File] filtered_summary_stats = summarize_filtered_sumstats.summary_stats_output
        Array[File] gzipped_filtered_summary_stats = summarize_filtered_sumstats.gzipped_summary_stats_output
        Array[Array[File]] sig_filtered_summary_stats = summarize_filtered_sumstats.sig_summary_stats_output
        Array[Array[Array[File]]] filtered_summary_stats_plots = summarize_filtered_sumstats.summary_plots
    }
}