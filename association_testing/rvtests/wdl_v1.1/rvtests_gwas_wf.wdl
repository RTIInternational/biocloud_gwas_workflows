version 1.1

import "rvtests_gwas_chr_wf.wdl" as RVCHR
import "collect_large_file_list_wf.wdl" as COLLECT
import "tsv_utils.wdl" as TSV
import "summarize_gwas_wf.wdl" as SUM
import "generate_kinship_matrix_wf.wdl" as KIN
import "convert_variant_ids.wdl" as IDCONVERT

workflow rvtests_gwas_wf{

    input{
        Array[File] vcfs_in
        Array[File] infos_in
        Array[String] chrs
        File pheno_file
        String pheno_name
        File? covar_file
        Array[String]? covars

        # Name of VCF annotation field that holds the allele dosage (usually 'DS')
        String dosage

        # Set output basename
        String study_output_basename

        # Boolean switches to select the correct pipeline depending on type of input data
        # Whether phenotype of interest is continous/binary
        Boolean is_continuous = false
        # Whether cohort includes unrelated or related individuals (and kinship matrix needs to be included in GWAS)
        Boolean is_related = false

        # Reference files for annotating output summary statistics wth 1000G population MAF info
        Array[File]? pop_maf_files
        String? maf_population

        # Imputation quality filtering cutoffs. Final outputs will have R2 >= this value
        Array[Float]? rsq_cutoffs

        # P-value cutoff for filtering set of interesting results
        Float sig_alpha

        # Optional MAF cutoffs levels used to generate summary stats/graphs across a range of MAF filters
        Array[Float]? maf_cutoffs

        # Boolean switches for whether to filter by sample/cohort MAF, 1000G populations MAF (can also be both/neither)
        Boolean filter_by_pop_maf = true
        Boolean filter_by_sample_maf = true

        # ID conversion parameters
        Array[File] id_ref_files
        String? in_missing_allele
        String? in_deletion_allele
        String ref_deletion_allele = "."
        Boolean rescue_rsids = false

        # RVTests arguments. These shouldn't change unless you REALLY know what you're doing.
        # Automatically set to the correct values based on continuous/binary flag above
        Boolean inverseNormal = is_continuous
        Boolean useResidualAsPhenotype = is_continuous
        Boolean qtl = is_continuous
        Array[String] metaTestsMaybe = ["score"]
        Boolean sex = false
        Boolean multipleAllele = false
        String? xLabel

        # Parameters for how/whether to further split VCFs for each chr
        # For extremely large VCFs this can speed up association testing
        # Cohorts under 5K individuals will probably run faster if you DONT split chr VCFs
        Boolean split_vcfs = false
        Int records_per_split = 100000
        Int split_vcf_cpu = 8

        # CPU/Memory each split of RVTests will consume
        # Keep this VERY SMALL (1CPU/2GB) if you're splitting within chr because RVTests is fast and you might have hundreds or thousands of splits
        Int rvtests_cpu_per_split = 1
        Int rvtests_mem_gb_per_split = 2

        # Options for kinship matrix (only used if is_realted == true)
        Boolean kinship_xHemi = true
        Boolean kinship_useBaldingNicols = true
        Boolean kinship_useIBS = false
        Float? kinship_minMAF = 0.05
        Float? kinship_minSiteQual
        Float? kinship_maxMiss
        # CPU/Mem to use for generating kinship matrix on each chr
        Int? split_kinship_cpu = 1
        Int? split_kinship_mem_gb = 2
        # CPU/Mem to use for the combineKinship task when generating a kinship matrix
        Int? combine_kinship_cpu = 2
        Int? combine_kinship_mem_gb = 4
        File? xHemiKinship
        File? kinship
        Boolean do_kinship_generate = (!defined(xHemiKinship)) || (!defined(kinship))

        # Mem for plotting
        Int plot_mem_gb

        # Container settings
        String image_source = "docker"
        String? ecr_repo
    }

    # Generate kinship matrix if samples are related
    if(is_related && do_kinship_generate){
        call KIN.generate_kinship_matrix_wf as get_kinship_mat{
            input:
                input_vcfs = vcfs_in,
                chrs = chrs,
                ped_file = pheno_file,
                output_basename = study_output_basename,
                xHemi = kinship_xHemi,
                useBaldingNicols = kinship_useBaldingNicols,
                useIBS = kinship_useIBS,
                dosage = dosage,
                xLabel = xLabel,
                maxMiss = kinship_maxMiss,
                minMAF = kinship_minMAF,
                minSiteQual = kinship_minSiteQual,
                split_kinship_cpu = split_kinship_cpu,
                split_kinship_mem_gb = split_kinship_mem_gb,
                combine_kinship_cpu = combine_kinship_cpu,
                combine_kinship_mem_gb = combine_kinship_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Choose kinship matrix to use
    # This is because you can either have no kinmatrix, have computed one in previous step, or passed one as input
    # Need additional logic because a simple select_first won't handle all cases
    File? kin_mat = if(defined(get_kinship_mat.kinship)) then get_kinship_mat.kinship else kinship
    File? xHemi_kin_mat = if(defined(get_kinship_mat.xHemiKinship)) then get_kinship_mat.xHemiKinship else xHemiKinship

    # Do RVTests chr workflow on each chromosome in parallel
    scatter(chr_index in range(length(vcfs_in))){

        pop_maf_file = select_first([pop_maf_files[chr_index], ""])
        call RVCHR.rvtests_gwas_chr_wf as rvtests{
            input:
                vcf_in = vcfs_in[chr_index],
                info_in = infos_in[chr_index],
                pheno_file = pheno_file,
                pheno_name = pheno_name,
                covar_file = covar_file,
                covars = covars,
                dosage = dosage,
                pop_maf_file = pop_maf_file,
                maf_population = maf_population,
                metaTestsMaybe = metaTestsMaybe,
                sex = sex,
                multipleAllele = multipleAllele,
                xLabel = xLabel,
                inverseNormal = inverseNormal,
                useResidualAsPhenotype = useResidualAsPhenotype,
                qtl = qtl,
                split_vcf_records = split_vcfs,
                records_per_split = records_per_split,
                split_vcf_cpu = split_vcf_cpu,
                rvtests_cpu = rvtests_cpu_per_split,
                rvtests_mem_gb = rvtests_mem_gb_per_split,
                kinship = kin_mat,
                xHemiKinship = xHemi_kin_mat,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Convert ids on summary stats output to standard ids
        call IDCONVERT.convert_variant_ids{
            input:
                chr = chrs[chr_index],
                in_file = rvtests.summary_stats,
                in_header = 1,
                in_sep = "tab",
                ref = id_ref_files[chr_index],
                in_id_col = 0,
                in_chr_col = 1,
                in_pos_col = 2,
                in_a1_col = 3,
                in_a2_col = 4,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = ref_deletion_allele,
                output_filename = basename(rvtests.summary_stats, ".tsv.gz") + ".good_ids.tsv.gz",
                rescue_rsids = rescue_rsids,
                output_compression = "gzip",
                cpu = rvtests_cpu_per_split,
                mem_gb = rvtests_mem_gb_per_split,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    # Merge chr outputs into single file
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
            input_files = convert_variant_ids.output_file,
            output_dir_name = study_output_basename + "_rvtests_chr_output",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = study_output_basename + ".rvtests.MetaAssoc.tsv",
            image_source = image_source,
            ecr_repo = ecr_repo
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
                sig_alpha = sig_alpha,
                sample_maf_cutoff = sample_maf_cutoff,
                pop_maf_cutoff = pop_maf_cutoff,
                min_rsq = min_rsq,
                plot_mem_gb = plot_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
            }
    }

    output{
        File raw_summary_stats = cat_sumstats.tsv_output
        Array[File] filtered_summary_stats = summarize_filtered_sumstats.summary_stats_output
        Array[File] gzipped_filtered_summary_stats = summarize_filtered_sumstats.gzipped_summary_stats_output
        Array[File] sig_filtered_summary_stats = summarize_filtered_sumstats.sig_summary_stats_output
        Array[Array[File]] filtered_summary_stats_plots = summarize_filtered_sumstats.summary_plots
    }
}