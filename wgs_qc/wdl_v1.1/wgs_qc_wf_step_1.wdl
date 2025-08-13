version 1.1

import "convert_variant_ids_bfile_chr_wf.wdl" as IDCONVERT
import "plink.wdl" as PLINK
import "remove_duplicate_variants_wf.wdl" as DUPLICATES
import "mahalanobis_ancestry_wf.wdl" as ANCESTRY
import "utils.wdl" as UTILS
import "tsv_utils.wdl" as TSV_UTILS
import "wgs_qc_wf_step_1_structs.wdl" as STRUCTS

workflow wgs_qc_wf_step_1{

    input {

        # Input dataset parameters
        Array[File] vcfs
        File keep
        String output_basename
        String dataset_short_name
        String dataset_display_name = dataset_short_name
        Array[String] chrs
        String genome_build_code

        # Ref files for ID conversion
        Array[File] id_conversion_ref_files

        # List of variants to use for ancestry and subject level QC
        File ref_variant_list

        # mahalanobis_ancestry_wf parameters
        File ancestry_ref_bed
        File ancestry_ref_bim
        File ancestry_ref_fam
        File ancestry_ref_psam
        String ancestry_pop_type = "SUPERPOP"
        Array[String] ancestries_to_include = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        Array[String] ancestries_display_names = ["African", "American Admixed", "East Asian", "European", "South Asian"]
        Int ancestry_std_dev_cutoff = 3

        # Sample filter cutoffs
        Float max_sample_missing_call_rate = 0.03
        Float min_sample_he = -0.2
        Float max_sample_he = 0.5

        # Variant filtering cutoffs
        Float max_variant_missing_call_rate = 0.03
        Float hwe_filter_pvalue = 0.0001

        # LD-Filtering parameters for subworkflows (smartpca_ancestry, relatedness, sex-check)
        File ld_exclude_regions
        String ld_type = "indep-pairwise"
        Int ld_window_size = 20000
        Int ld_step_size = 2000
        Float ld_r2_threshold = 0.5
        Float ld_maf_cutoff = 0.01
        
        # Cutoff below which an ancestry group of samples won't go through full pipeline
        # Handles cases where you might only be excluding a handful of outlier samples and only care about the main ancestry groups
        Int min_ancestry_samples_to_postprocess = 50

        # Sex-check filter parameters
        Boolean filter_discrepant_sex = true
        Float max_female_f = 0.2
        Float min_male_f = 0.8

        # Sex xref for adding sex
        File sex_xref
        Int sex_xref_fid_col     # 1-based
        Int sex_xref_iid_col     # 1-based
        Int sex_xref_sex_col     # 1-based
        Boolean sex_xref_header = false
        String? sex_xref_delimiter

        # Kinship filter parameters
        Boolean filter_related_samples = true
        Int degree = 3
        Int num_king_splits = 4

        # PCA parameters
        Int kinship_pcs_to_analyze = 3

        # Ancestral SNP filtering parameters
        Float ancestral_pca_loading_cutoff = 0.003
        Int max_kinship_snps = 100000
        Int min_kinship_snps = 10000
        Float ancestral_pca_loading_step_size = 0.001
        Float max_ancestral_pca_loading_cutoff = 0.01

        # Container
        String image_source = "docker"
        String? ecr_repo
    
    }

    # Check some common errors to save time because otherwise it would take like 2 hours to catch these
    # Quit if ancestry definitions and ancestries aren't same length
    if(length(ancestries_to_include) != length(ancestries_display_names)){
        call UTILS.raise_error as ancestry_definition_fail{
            input:
                msg = "Workflow input error: ancestries_to_incude must be same length as ancestries_display_names!",
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Quit if ancestry != POP or SUPERPOP
    if(ancestry_pop_type != "SUPERPOP" && ancestry_pop_type != "POP"){
        call UTILS.raise_error as ancestry_pop_type_fail{
            input:
                msg = "Workflow input error: ancestry_type MUST be either POP or SUPERPOP!",
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Prepare file for updating sex
    Array[String] sex_xref_cols = ["~{sex_xref_fid_col}", "~{sex_xref_iid_col}", "~{sex_xref_sex_col}"]
    call TSV_UTILS.tsv_select as create_sex_update_file {
        input:
            tsv_input = sex_xref,
            output_filename = "sex_update.tsv",
            fields = sex_xref_cols,
            header = sex_xref_header,
            delimiter = sex_xref_delimiter,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    scatter(i in range(length(chrs))) {

        String chr = "~{chrs[i]}"

        Int vcf_multiplier = floor((size(vcfs[i], "GB") / 5) + 1)
        
        # Convert vcf to plink bfile
        call PLINK.convert_vcf_to_bed as all_snps_chr_convert_vcf_to_bfile{
            input:
                vcf_in = vcfs[i],
                output_basename = "~{output_basename}_chr~{chr}",
                cpu = 4 * vcf_multiplier,
                mem_gb = 8 * vcf_multiplier,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Get initial sample count
        call UTILS.wc as sample_init_count {
            input:
                input_file = keep,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
        Int sample_multiplier_chr = floor((sample_init_count.num_lines / 750) + 1)

        # Get initial variant count
        call UTILS.wc as variant_init_count {
            input:
                input_file = all_snps_chr_convert_vcf_to_bfile.bim_out,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
        Int variant_multiplier_chr = floor((variant_init_count.num_lines / 7500000) + 1)

        Int combined_multiplier_chr = sample_multiplier_chr * variant_multiplier_chr

        # Extract batch samples and update sex
        call PLINK.make_bed as all_snps_chr_select_samples{
            input:
                bed_in = all_snps_chr_convert_vcf_to_bfile.bed_out,
                bim_in = all_snps_chr_convert_vcf_to_bfile.bim_out,
                fam_in = all_snps_chr_convert_vcf_to_bfile.fam_out,
                keep_samples = keep,
                update_sex = create_sex_update_file.tsv_output,
                output_basename = "~{output_basename}_chr~{chr}_batch_samples",
                cpu = floor((0.5 * combined_multiplier_chr) + 1),
                mem_gb = 2 * combined_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Remove phenotype from fam file
        call PLINK.remove_fam_phenotype as remove_fam_pheno{
            input:
                fam_in = all_snps_chr_select_samples.fam_out,
                output_basename = "~{output_basename}_chr~{chr}_no_pheno",
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Fix founder information in pedigree to make sure mother and father ids actually exist in dataset
        call PLINK.make_founders as make_founders{
            input:
                fam_in = remove_fam_pheno.fam_out,
                output_basename = "~{output_basename}_chr~{chr}_no_pheno_make_founders",
                mem_gb = 1 * sample_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Convert variant IDs to standard format
        call IDCONVERT.convert_variant_ids_bfile_chr_wf as all_snps_chr_convert_variant_ids{
            input:
                bim_in = all_snps_chr_select_samples.bim_out,
                ref_file = id_conversion_ref_files[i],
                chr = chr,
                output_basename = "~{output_basename}_chr~{chr}_convert_ids.bim",
                in_chunk_size = 50000,
                ref_chunk_size = 1000000,
                convert_cpu = 2,
                convert_mem_gb = 4,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Remove duplicates
        call DUPLICATES.remove_duplicate_variants_wf as all_snps_chr_remove_duplicates{
            input:
                bed_in = all_snps_chr_select_samples.bed_out,
                bim_in = all_snps_chr_convert_variant_ids.bim_out,
                fam_in = make_founders.fam_out,
                output_basename = "~{output_basename}_chr~{chr}_no_dups",
                plink_mem_gb = 3 * combined_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Extract ref_panel variants from dataset
        call PLINK.make_bed as ref_snps_chr_extract_variants{
            input:
                bed_in = all_snps_chr_remove_duplicates.bed_out,
                bim_in = all_snps_chr_remove_duplicates.bim_out,
                fam_in = all_snps_chr_remove_duplicates.fam_out,
                extract = ref_variant_list,
                output_basename = "~{output_basename}_chr~{chr}_ref_snps",
                cpu = floor((0.25 * combined_multiplier_chr) + 1),
                mem_gb = floor((0.5 * combined_multiplier_chr) + 1),
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    # Get multiplier for whole genome
    call UTILS.sum_ints as sum_variant_multipliers{
        input:
            ints = variant_multiplier_chr,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int combined_multiplier = sample_multiplier_chr[0] * sum_variant_multipliers.sum

    # Merge chromosomes for ref panel variants
    call PLINK.merge_beds as ref_snps_post_id_conversion{
        input:
            bed_in = ref_snps_chr_extract_variants.bed_out,
            bim_in = ref_snps_chr_extract_variants.bim_out,
            fam_in = ref_snps_chr_extract_variants.fam_out,
            output_basename = "~{output_basename}_ref_snps",
            cpu = floor((combined_multiplier / 50) + 1),
            mem_gb = floor((combined_multiplier / 2) + 1),
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Smartpca ancestry WF to partition by ancestry
    call ANCESTRY.mahalanobis_ancestry_wf{
        input:
            dataset_bed = ref_snps_post_id_conversion.bed_out,
            dataset_bim = ref_snps_post_id_conversion.bim_out,
            dataset_fam = ref_snps_post_id_conversion.fam_out,
            dataset_short_name = dataset_short_name,
            dataset_display_name = dataset_display_name,
            ref_bed = ancestry_ref_bed,
            ref_bim = ancestry_ref_bim,
            ref_fam = ancestry_ref_fam,
            ref_psam = ancestry_ref_psam,
            ancestry_pop_type = ancestry_pop_type,
            ancestries_to_include = ancestries_to_include,
            ancestries_display_names = ancestries_display_names,
            do_id_conversion = false,
            ld_exclude_regions = ld_exclude_regions,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Merge chr bims
    Array[Array[File]] bims_to_merge = [
        all_snps_chr_convert_variant_ids.bim_out,
        all_snps_chr_remove_duplicates.bim_out
    ]
    Array[String] merged_bim_filenames = [
        "~{output_basename}_init.bim",
        "~{output_basename}_no_dups.bim"
    ]
    scatter(i in range(length(bims_to_merge))) {
        call UTILS.cat as merge_bims{
            input:
                input_files = bims_to_merge[i],
                output_filename = merged_bim_filenames[i],
                mem_gb = floor((sum_variant_multipliers.sum / 5) + 1),
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get SNP lists
    Array[File] bims_for_snp_lists = [
        merge_bims.output_file[0],
        merge_bims.output_file[1]
    ]
    Array[String] snp_list_filenames = [
        "~{output_basename}_variants_initial.txt",
        "~{output_basename}_variants_final.txt",
    ]
    scatter(i in range(length(bims_for_snp_lists))) {
        call UTILS.cut as snp_list{
            input:
                input_file = bims_for_snp_lists[i],
                args =  "-f 2",
                output_filename = snp_list_filenames[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get list of duplicate SNPs
    call UTILS.comm as duplicate_snps{
        input:
            file1 = snp_list.output_file[0],
            file2 = snp_list.output_file[1],
            option = "-23",
            output_filename = "~{output_basename}_variants_duplicates.txt",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get SNP counts
    Array[File] files_for_snp_counts = [
        snp_list.output_file[0],
        snp_list.output_file[1],
        duplicate_snps.output_file
    ]
    scatter(i in range(length(files_for_snp_counts))) {
        call UTILS.wc as snp_count {
            input:
                input_file = files_for_snp_counts[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get sample lists
    Array[File] fams_for_sample_lists = [
        all_snps_chr_convert_vcf_to_bfile.fam_out[0],
        all_snps_chr_select_samples.fam_out[0]
    ]
    Array[String] sample_list_filenames = [
        "~{output_basename}_samples_initial.tsv",
        "~{output_basename}_samples_final.tsv"
    ]
    scatter(i in range(length(fams_for_sample_lists))) {
        call UTILS.cut as sample_list{
            input:
                input_file = fams_for_sample_lists[i],
                args =  "-f 1,2",
                output_filename = sample_list_filenames[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }
    Array[Pair[String, File]] ancestry_samples = zip(mahalanobis_ancestry_wf.ancestry_labels, mahalanobis_ancestry_wf.dataset_ancestry_keep_lists)

    # Get sample counts
    Array[File] files_for_sample_counts = [
        sample_list.output_file[0],
        sample_list.output_file[1]
    ]
    scatter(i in range(length(files_for_sample_counts))) {
        call UTILS.wc as sample_count {
            input:
                input_file = files_for_sample_counts[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }
    scatter(i in range(length(mahalanobis_ancestry_wf.dataset_ancestry_keep_lists))) {
        call UTILS.wc as ancestry_sample_count {
            input:
                input_file = mahalanobis_ancestry_wf.dataset_ancestry_keep_lists[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }
    Array[Pair[String, Int]] ancestry_sample_counts = zip(mahalanobis_ancestry_wf.ancestry_labels, ancestry_sample_count.num_lines)

    output{

        # Genotypes
        GENOTYPE_FILES genotype_files = GENOTYPE_FILES {
            beds: all_snps_chr_remove_duplicates.bed_out,
            bims: all_snps_chr_remove_duplicates.bim_out,
            fams: all_snps_chr_remove_duplicates.fam_out
        }

        # Pruned genotypes
        PRUNED_GENOTYPE_FILES pruned_genotype_files = PRUNED_GENOTYPE_FILES {
            bed: ref_snps_post_id_conversion.bed_out,
            bim: ref_snps_post_id_conversion.bim_out,
            fam: ref_snps_post_id_conversion.fam_out
        }

        # Variant lists
        VARIANT_LISTS variant_lists = VARIANT_LISTS {
            initial: snp_list.output_file[0],
            final: snp_list.output_file[1],
            duplicates: duplicate_snps.output_file
        }

        # Sample lists
        SAMPLE_LISTS sample_lists = SAMPLE_LISTS {
            initial: sample_list.output_file[0],
            final: sample_list.output_file[1],
            ancestries: ancestry_samples
        }

        # Counts
        COUNTS counts = COUNTS {
            variant_initial: snp_count.num_lines[0],
            variant_final: snp_count.num_lines[1],
            variant_duplicates: snp_count.num_lines[2],
            sample_initial: sample_count.num_lines[0],
            sample_final: sample_count.num_lines[1],
            sample_ancestries: ancestry_sample_counts
        }

        # Outputs from ancestry assignment workflow
        ANCESTRY_WF_OUTPUTS ancestry_wf_outputs = ANCESTRY_WF_OUTPUTS {
            dataset_ancestry_labels: mahalanobis_ancestry_wf.ancestry_labels,
            dataset_ancestry_keep_lists: mahalanobis_ancestry_wf.dataset_ancestry_keep_lists,
            dataset_ancestry_assignments: mahalanobis_ancestry_wf.dataset_ancestry_assignments,
            dataset_ancestry_assignments_summary: mahalanobis_ancestry_wf.dataset_ancestry_assignments_summary,
            dataset_ancestry_assignments_plots: mahalanobis_ancestry_wf.dataset_ancestry_assignments_plots,
            evec: mahalanobis_ancestry_wf.evec,
            eval: mahalanobis_ancestry_wf.eval,
            snpweight: mahalanobis_ancestry_wf.snpweight,
            log: mahalanobis_ancestry_wf.smartpca_log,
            ref_dropped_samples: mahalanobis_ancestry_wf.ref_dropped_samples,
            ref_raw_ancestry_assignments: mahalanobis_ancestry_wf.ref_raw_ancestry_assignments,
            ref_raw_ancestry_assignments_summary: mahalanobis_ancestry_wf.ref_raw_ancestry_assignments_summary,
            pre_processing_pc_plots: mahalanobis_ancestry_wf.pre_processing_pc_plots,
            dataset_ancestry_outliers_plots: mahalanobis_ancestry_wf.dataset_ancestry_outliers_plots
        }

        # Step 2 parameters
        STEP_2_PARAMETERS step_2_parameters = STEP_2_PARAMETERS {
            complete_beds: all_snps_chr_remove_duplicates.bed_out,
            complete_bims: all_snps_chr_remove_duplicates.bim_out,
            complete_fams: all_snps_chr_remove_duplicates.fam_out,
            filtered_bed: ref_snps_post_id_conversion.bed_out,
            filtered_bim: ref_snps_post_id_conversion.bim_out,
            filtered_fam: ref_snps_post_id_conversion.fam_out,
            output_basename: output_basename,
            chrs: chrs,
            genome_build_code: genome_build_code,
            variant_init_count: snp_count.num_lines[1],
            sample_init_count: sample_count.num_lines[1],
            ld_exclude_regions: ld_exclude_regions,
            ld_type: ld_type,
            ld_window_size: ld_window_size,
            ld_step_size: ld_step_size,
            ld_r2_threshold: ld_r2_threshold,
            ld_maf_cutoff: ld_maf_cutoff,
            max_variant_missing_call_rate: max_variant_missing_call_rate,
            hwe_filter_pvalue: hwe_filter_pvalue,
            max_sample_missing_call_rate: max_sample_missing_call_rate,
            min_sample_he: min_sample_he,
            max_sample_he: max_sample_he,
            filter_related_samples: filter_related_samples,
            degree: degree,
            num_king_splits: num_king_splits,
            kinship_pcs_to_analyze: kinship_pcs_to_analyze,
            filter_discrepant_sex: filter_discrepant_sex,
            max_female_f: max_female_f,
            min_male_f: min_male_f,
            ancestral_pca_loading_cutoff: ancestral_pca_loading_cutoff,
            max_kinship_snps: max_kinship_snps,
            min_kinship_snps: min_kinship_snps,
            ancestral_pca_loading_step_size: ancestral_pca_loading_step_size,
            max_ancestral_pca_loading_cutoff: max_ancestral_pca_loading_cutoff,
            image_source: image_source,
            ecr_repo: ecr_repo
        }

    }

}
