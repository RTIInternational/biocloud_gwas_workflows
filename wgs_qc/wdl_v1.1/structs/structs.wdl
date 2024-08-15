version 1.1

struct STEP_2_ARGS {
    String ancestry
    Array[File] complete_beds
    Array[File] complete_bims
    Array[File] complete_fams
    File filtered_bed
    File filtered_bim
    File filtered_fam
    File ancestry_samples
    String output_basename
    Array[String] chrs
    String genome_build_code
    Float max_missing_site_rate
    Float hwe_filter_pvalue
    File? ld_exclude_regions
    String ld_type
    Int ld_window_size
    Int ld_step_size
    Float ld_r2_threshold
    Float ld_maf_cutoff
    Float max_sample_missing_rate
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
    Int plink_cpu
    Int plink_mem_gb
    Int plink_chr_cpu
    Int plink_chr_mem_gb
    Int merge_bed_cpu
    Int merge_bed_mem_gb
    Int sex_check_cpu
    Int sex_check_mem_gb
    Int king_cpu_per_split
    Int king_mem_gb_per_split
    Int pca_cpu
    Int pca_mem_gb
    String image_source
    String ecr_repo
}

struct COUNTS {
    Int variant_count_initial
    Int variant_count_duplicates
    Int variant_count_no_duplicates
    Int sample_count_initial
    Int sample_count_batch
    Array[Pair[String, Int]] sample_count_ancestry
}
