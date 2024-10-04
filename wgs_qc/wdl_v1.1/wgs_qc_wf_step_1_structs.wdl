version 1.1

struct GENOTYPE_FILES {
    Array[File] beds
    Array[File] bims
    Array[File] fams
}

struct PRUNED_GENOTYPE_FILES {
    File bed
    File bim
    File fam
}

struct VARIANT_LISTS {
    File initial
    File final
    File duplicates
}

struct SAMPLE_LISTS {
    File initial
    File final
    Array[Pair[String, File]] ancestries
}

struct COUNTS {
    Int variant_initial
    Int variant_final
    Int variant_duplicates
    Int sample_initial
    Int sample_final
    Array[Pair[String, Int]] sample_ancestries
}

struct ANCESTRY_WF_OUTPUTS {
    Array[String] dataset_ancestry_labels
    Array[File] dataset_ancestry_keep_lists
    File dataset_ancestry_assignments
    File dataset_ancestry_assignments_summary
    Array[File] dataset_ancestry_assignments_plots
    File evec
    File eval
    File snpweight
    File log
    File ref_dropped_samples
    File ref_raw_ancestry_assignments
    File ref_raw_ancestry_assignments_summary
    Array[File] pre_processing_pc_plots
    Array[File] dataset_ancestry_outliers_plots
}

struct STEP_2_PARAMETERS {
    Array[File] complete_beds
    Array[File] complete_bims
    Array[File] complete_fams
    File filtered_bed
    File filtered_bim
    File filtered_fam
    String output_basename
    Array[String] chrs
    String genome_build_code
    Int variant_init_count
    Int sample_init_count
    File ld_exclude_regions
    String ld_type
    Int ld_window_size
    Int ld_step_size
    Float ld_r2_threshold
    Float ld_maf_cutoff
    Float max_variant_missing_call_rate
    Float hwe_filter_pvalue
    Float max_sample_missing_call_rate
    Float min_sample_he
    Float max_sample_he
    Boolean filter_related_samples
    Int degree
    Int num_king_splits
    Int kinship_pcs_to_analyze
    Boolean filter_discrepant_sex
    Float max_female_f
    Float min_male_f
    Float ancestral_pca_loading_cutoff
    Int max_kinship_snps
    Int min_kinship_snps
    Float ancestral_pca_loading_step_size
    Float max_ancestral_pca_loading_cutoff
    String image_source
    String ecr_repo
}
