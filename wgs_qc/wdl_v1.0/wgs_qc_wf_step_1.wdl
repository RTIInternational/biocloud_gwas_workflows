version 1.0

import "convert_variant_ids_bfile_chr_wf.wdl" as IDCONVERT
import "plink.wdl" as PLINK
import "remove_duplicate_variants_wf.wdl" as DUPLICATES
import "smartpca_ancestry_wf.wdl" as ANCESTRY
import "utils.wdl" as UTILS
import "tsv_utils.wdl" as TSV_UTILS
import "structs.wdl" as STRUCTS

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

        #### Parameters for id conversion
        Array[File] id_conversion_ref_files
        String? id_conversion_in_sep
        String? id_conversion_in_missing_allele
        String? id_conversion_in_deletion_allele
        String? id_conversion_ref_deletion_allele
        Int? id_conversion_in_chunk_size
        Int? id_conversion_ref_chunk_size
        Boolean? id_conversion_rescue_rsids

        # List of variants to use for ancestry and subject level QC
        File ref_variant_list

        # smartpca_ancestry_wf parameters
        File ancestry_ref_bed
        File ancestry_ref_bim
        File ancestry_ref_fam
        File ancestry_ref_psam
        String ancestry_pop_type = "SUPERPOP"
        Array[String] ancestries_to_include = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        Array[String] ancestries_display_names = ["African", "American Admixed", "East Asian", "European", "South Asian"]
        Int ancestry_std_dev_cutoff = 3

        # Sample filter cutoffs
        Float max_sample_missing_rate = 0.03
        Float min_sample_he = -0.2
        Float max_sample_he = 0.5

        # Variant filtering cutoffs
        Float max_missing_site_rate = 0.03
        Float hwe_filter_pvalue = 0.0001

        # LD-Filtering parameters for subworkflows (smartpca_ancestry, relatedness, sex-check)
        File? ld_exclude_regions
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

        # Resources for various modules

        # General CPU/Mem for basically any PLINK task that runs on the whole dataset
        Int plink_cpu = 1
        Int plink_mem_gb = 2

        # General CPU/MEM for basically any PLINK task that runs on a single chr
        Int plink_chr_cpu = 1
        Int plink_chr_mem_gb = 2

        # General CPU/MEM for jobs that merge bed files either by chr or sample
        Int merge_bed_cpu = 4
        Int merge_bed_mem_gb = 8

        # CPU/MEM for conversion from vcf to bfile
        Int plink_vcf_to_bed_cpu = 8
        Int plink_vcf_to_bed_mem_gb = 64

        # Resource for converting ids on each chr
        Int id_conversion_cpu = 1
        Int id_conversion_mem_gb = 4

        # Speicific tasks where resource limits may need to be adjusted for larger/smaller inputs
        Int sex_check_cpu = 4
        Int sex_check_mem_gb = 8
        Int king_cpu_per_split = 4
        Int king_mem_gb_per_split = 8
        Int pca_cpu = 4
        Int pca_mem_gb = 8
        Int assign_ancestry_cpu = 4
        Int assign_ancestry_mem_gb = 8

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

        # Convert vcf to plink bfile
        call PLINK.convert_vcf_to_bed as all_snps_chr_convert_vcf_to_bfile{
            input:
                vcf_in = vcfs[i],
                output_basename = "~{output_basename}_chr~{chr}",
                cpu = plink_vcf_to_bed_cpu,
                mem_gb = plink_vcf_to_bed_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        call UTILS.wc as chr_count_init_snps{
            input:
                input_file = all_snps_chr_convert_vcf_to_bfile.bim_out,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Extract batch samples and update sex
        call PLINK.make_bed as all_snps_chr_select_samples{
            input:
                bed_in = all_snps_chr_convert_vcf_to_bfile.bed_out,
                bim_in = all_snps_chr_convert_vcf_to_bfile.bim_out,
                fam_in = all_snps_chr_convert_vcf_to_bfile.fam_out,
                keep_samples = keep,
                update_sex = create_sex_update_file.tsv_output,
                output_basename = "~{output_basename}_chr~{chr}_batch_samples",
                cpu = plink_chr_cpu,
                mem_gb = plink_chr_mem_gb,
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
                cpu = plink_chr_cpu,
                mem_gb = plink_chr_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }


        # Convert variant IDs to standard format
        call IDCONVERT.convert_variant_ids_bfile_chr_wf as all_snps_chr_convert_variant_ids{
            input:
                bim_in = all_snps_chr_select_samples.bim_out,
                ref_file = id_conversion_ref_files[i],
                chr = chr,
                output_filename = "~{output_basename}_chr~{chr}_convert_ids.bim",
                in_sep = id_conversion_in_sep,
                in_missing_allele = id_conversion_in_missing_allele,
                in_deletion_allele = id_conversion_in_deletion_allele,
                ref_deletion_allele = id_conversion_ref_deletion_allele,
                in_chunk_size = id_conversion_in_chunk_size,
                ref_chunk_size = id_conversion_ref_chunk_size,
                rescue_rsids = id_conversion_rescue_rsids,
                convert_cpu = id_conversion_cpu,
                convert_mem_gb = id_conversion_mem_gb,
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
                label_duplicate_variants_cpu = id_conversion_cpu,
                label_duplicate_variants_mem_gb = id_conversion_mem_gb,
                plink_cpu = plink_chr_cpu,
                plink_mem_gb = plink_chr_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        call UTILS.wc as chr_count_post_id_conversion_snps{
            input:
                input_file = all_snps_chr_remove_duplicates.bim_out,
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
                cpu = plink_chr_cpu,
                mem_gb = plink_chr_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    # Count initial snps in dataset
    call UTILS.sum_ints as count_init_snps{
        input: 
            ints = chr_count_init_snps.num_lines,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int init_snp_count = count_init_snps.sum

    # Count initial samples in batch
    call UTILS.wc as count_init_samples{
        input:
            input_file = all_snps_chr_select_samples.fam_out[0],
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int init_sample_count = count_init_samples.num_lines

    # Count post ID conversion snps
    call UTILS.sum_ints as count_post_id_conversion_snps{
        input: 
            ints = chr_count_post_id_conversion_snps.num_lines,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int post_id_conversion_snp_count = count_post_id_conversion_snps.sum

    # Merge chromosomes for ref panel variants
    call PLINK.merge_beds as ref_snps_post_id_conversion{
        input:
            bed_in = ref_snps_chr_extract_variants.bed_out,
            bim_in = ref_snps_chr_extract_variants.bim_out,
            fam_in = ref_snps_chr_extract_variants.fam_out,
            output_basename = "~{output_basename}_ref_snps",
            cpu = merge_bed_cpu,
            mem_gb = merge_bed_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Smartpca ancestry WF to partition by ancestry
    call ANCESTRY.smartpca_ancestry_wf{
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
            id_conversion_ref_files = id_conversion_ref_files,
            ld_exclude_regions = ld_exclude_regions,
            plink_dataset_cpu = plink_cpu,
            plink_dataset_mem_gb = plink_mem_gb,
            convert_variant_ids_cpu = id_conversion_cpu,
            convert_variant_ids_mem_gb = id_conversion_mem_gb,
            smartpca_cpu = pca_cpu,
            smartpca_mem_gb = pca_mem_gb,
            assign_ancestry_cpu = assign_ancestry_cpu,
            assign_ancestry_mem_gb = assign_ancestry_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Filter out any ancestries with less than a minimum cutoff of samples
    scatter(ancestry_index in range(length(ancestries_to_include))){
        # Count number of samples in each ancestry
        call UTILS.wc as count_ancestry_samples{
            input:
                input_file = smartpca_ancestry_wf.dataset_ancestry_keep_lists[ancestry_index],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
        Int ancestry_sample_count = count_ancestry_samples.num_lines

        # Only include ancestries exceeding minimum number of samples
        if(ancestry_sample_count > min_ancestry_samples_to_postprocess){
            Int ancestry_sample_count_maybe = ancestry_sample_count
            String ancestry_maybe = ancestries_to_include[ancestry_index]
            File ancestry_samples_maybe = smartpca_ancestry_wf.dataset_ancestry_keep_lists[ancestry_index]
        }
    }

    # Remove null values from filtered ancestries to create list of ancestries to postprocess
    Array[Int] input_ancestries_init_sample_counts = ancestry_sample_count
    Array[Int] out_ancestries_init_sample_counts = select_all(ancestry_sample_count_maybe)
    Array[String] out_ancestries = select_all(ancestry_maybe)
    Array[File] output_ancestry_samples = select_all(ancestry_samples_maybe)

    # Create JSON files with arguments for step 2 of WF
    scatter(ancestry_index in range(length(out_ancestries))){

        String ancestry = out_ancestries[ancestry_index]
        File ancestry_samples = output_ancestry_samples[ancestry_index]
        Int initial_ancestry_sample_count = out_ancestries_init_sample_counts[ancestry_index]
        STEP_2_ARGS step_2_args = object {
            ancestry: ancestry,
            complete_beds: all_snps_chr_remove_duplicates.bed_out,
            complete_bims: all_snps_chr_remove_duplicates.bim_out,
            complete_fams: all_snps_chr_remove_duplicates.fam_out,
            filtered_bed: ref_snps_post_id_conversion.bed_out,
            filtered_bim: ref_snps_post_id_conversion.bim_out,
            filtered_fam: ref_snps_post_id_conversion.fam_out,
            ancestry_samples: ancestry_samples,
            output_basename: "~{output_basename}_~{ancestry}",
            chrs: chrs,
            genome_build_code: genome_build_code,
            init_snp_count: init_snp_count,
            post_id_conversion_snp_count: post_id_conversion_snp_count,
            initial_ancestry_sample_count: initial_ancestry_sample_count,
            max_missing_site_rate: max_missing_site_rate,
            hwe_filter_pvalue: hwe_filter_pvalue,
            ld_exclude_regions: ld_exclude_regions,
            ld_type: ld_type,
            ld_window_size: ld_window_size,
            ld_step_size: ld_step_size,
            ld_r2_threshold: ld_r2_threshold,
            ld_maf_cutoff: ld_maf_cutoff,
            max_sample_missing_rate: max_sample_missing_rate,
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
            plink_cpu: plink_cpu,
            plink_mem_gb: plink_mem_gb,
            plink_chr_cpu: plink_chr_cpu,
            plink_chr_mem_gb: plink_chr_mem_gb,
            merge_bed_cpu: merge_bed_cpu,
            merge_bed_mem_gb: merge_bed_mem_gb,
            sex_check_cpu: sex_check_cpu,
            sex_check_mem_gb: sex_check_mem_gb,
            king_cpu_per_split: king_cpu_per_split,
            king_mem_gb_per_split: king_mem_gb_per_split,
            pca_cpu: pca_cpu,
            pca_mem_gb: pca_mem_gb,
            image_source: image_source,
            ecr_repo: ecr_repo
        }


        call write_step_2_args_to_json {
            input:
                json_file = "~{output_basename}_~{ancestry}_step_2_args.json",
                step_2_args = step_2_args,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    output{

        # Step 2 args
        Array[File] step_2_args_jsons = write_step_2_args_to_json.step_2_args_json

        # Filter count metrics
        Int initial_snp_count = init_snp_count
        Int initial_sample_count = init_sample_count
        Int duplicate_snp_count = init_snp_count - post_id_conversion_snp_count
        Array[String] input_ancestries = ancestries_to_include
        Array[Int] input_ancestries_initial_sample_counts = input_ancestries_init_sample_counts
        Array[String] output_ancestries = out_ancestries
        Array[Int] output_ancestries_initial_sample_counts = out_ancestries_init_sample_counts

        # Genotypes
        Array[File] full_beds = all_snps_chr_remove_duplicates.bed_out
        Array[File] full_bims = all_snps_chr_remove_duplicates.bim_out
        Array[File] full_fams = all_snps_chr_remove_duplicates.fam_out
        File pruned_bed = ref_snps_post_id_conversion.bed_out
        File pruned_bim = ref_snps_post_id_conversion.bim_out
        File pruned_fam = ref_snps_post_id_conversion.fam_out

        # Outputs from ancestry assignment workflow
        File ancestry_wf_evec = smartpca_ancestry_wf.evec
        File ancestry_wf_eval = smartpca_ancestry_wf.eval
        File ancestry_wf_snpweight = smartpca_ancestry_wf.snpweight
        File ancestry_wf_log = smartpca_ancestry_wf.smartpca_log
        Array[File] ancestry_wf_pre_processing_pc_plots  = smartpca_ancestry_wf.pre_processing_pc_plots
        File ancestry_wf_ref_dropped_samples  = smartpca_ancestry_wf.ref_dropped_samples
        File ancestry_wf_ref_raw_ancestry_assignments = smartpca_ancestry_wf.ref_raw_ancestry_assignments
        File ancestry_wf_ref_raw_ancestry_assignments_summary = smartpca_ancestry_wf.ref_raw_ancestry_assignments_summary
        File ancestry_wf_dataset_ancestry_assignments = smartpca_ancestry_wf.dataset_ancestry_assignments
        File ancestry_wf_dataset_ancestry_assignments_summary = smartpca_ancestry_wf.dataset_ancestry_assignments_summary
        Array[File] ancestry_wf_dataset_ancestry_assignments_plots = smartpca_ancestry_wf.dataset_ancestry_assignments_plots
        Array[File] ancestry_wf_dataset_ancestry_outliers_plots = smartpca_ancestry_wf.dataset_ancestry_outliers_plots
        Array[File] ancestry_wf_dataset_ancestry_keep_lists = smartpca_ancestry_wf.dataset_ancestry_keep_lists

    }
}


task write_step_2_args_to_json{

    input{

        String json_file
        STEP_2_ARGS step_2_args
        
        # Runtime environment
        String docker_image = "python:3.11.9-slim-bookworm"
        String ecr_image = "rtibiocloud/python:3.11.9-slim-bookworm"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 1

    }

    command <<<
        python <<CODE

        import json
        
        with open("~{write_json(test_struct_instance)}") as f:
            test_struct = json.load(f)
        
        with open("~{json_file}", "w") as f:
            json.dump(test_struct, f)
        
        CODE
    >>>


    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File step_2_args_json = "~{json_file}"
    }

}
