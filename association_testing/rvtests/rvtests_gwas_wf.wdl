import "biocloud_gwas_workflows/association_testing/rvtests/rvtests_gwas_chr_wf.wdl" as RVCHR
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/summarize_gwas_wf.wdl" as SUM
import "biocloud_gwas_workflows/helper_workflows/generate_kinship_matrix_wf.wdl" as KIN

workflow rvtests_gwas_wf{
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
    Boolean is_continuous
    # Whether cohort includes unrelated or related individuals (and kinship matrix needs to be included in GWAS)
    Boolean is_related

    # Reference files for annotating output summary statistics wth 1000G population MAF info
    Array[File] pop_maf_files
    String maf_population

    # Imputation quality filtering cutoff. Final outputs will have R2 >= this value
    Float min_rsq

    # P-value cutoff for filtering set of interesting results
    Float sig_alpha

    # Optional MAF cutoffs levels used to generate summary stats/graphs across a range of MAF filters
    Array[Float]? maf_cutoffs
    # Boolean switches for whether to filter by sample/cohort MAF, 1000G populations MAF (can also be both/neither)
    Boolean filter_by_pop_maf
    Boolean filter_by_sample_maf


    # RVTests arguments. These shouldn't change unless you REALLY know what you're doing.
    # Automatically set to the correct values based on continuous/binary flag above
    Boolean inverseNormal = is_continuous
    Boolean useResidualAsPhenotype = is_continuous
    Boolean qtl = is_continuous
    Array[String] metaTestsMaybe = ["score"]
    Boolean? sex
    Boolean? multipleAllele
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
    Boolean? kinship_xHemi = true
    Boolean? kinship_useBaldingNicols = true
    Boolean? kinship_useIBS = false
    Float? kinship_minMAF = 0.05
    Float? kinship_minSiteQual
    Float? kinship_maxMiss
    # CPU/Mem to use for generating kinship matrix on each chr
    Int? split_kinship_cpu = 1
    Int? split_kinship_mem_gb = 2
    # CPU/Mem to use for the combineKinship task when generating a kinship matrix
    Int? combine_kinship_cpu = 2
    Int? combine_kinship_mem_gb = 4

    # Generate kinship matrix if samples are related
    if(is_related){
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
                combine_kinship_mem_gb = combine_kinship_mem_gb
        }
    }

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
                kinship = get_kinship_mat.kinship,
                xHemiKinship = get_kinship_mat.xHemiKinship
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
                filter_string = "--is-numeric '10' --ge '10:${min_rsq}'"
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