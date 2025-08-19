version 1.1

import "plink.wdl" as PLINK
import "hwe_filter_chr_wf.wdl" as HWE
import "relatedness_wf.wdl" as REL
import "sex_check_wf.wdl" as SEX
import "utils.wdl" as UTILS
import "wgs_qc_wf_step_2_structs.wdl" as STRUCTS

workflow wgs_qc_ancestry_wf_step_2{

    input{

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

        # Initial counts
        Int variant_init_count
        Int sample_init_count

        # Variant filtering cutoffs
        Float max_variant_missing_call_rate = 0.03
        Float hwe_filter_pvalue = 0.0001

        # LD-Filtering parameters for subworkflows (smartpca_ancestry, relatedness, sex-check)
        File? ld_exclude_regions
        String ld_type = "indep-pairwise"
        Int ld_window_size = 20000
        Int ld_step_size = 2000
        Float ld_r2_threshold = 0.5
        Float ld_maf_cutoff = 0.01

        # Sample filter cutoffs
        Float max_sample_missing_call_rate = 0.03
        Float min_sample_he = -0.2
        Float max_sample_he = 0.5

        # Kinship filter parameters
        Boolean filter_related_samples = true
        Int degree = 3
        Int num_king_splits = 4

        # PCA parameters
        Int kinship_pcs_to_analyze = 3

        # Sex-check filter parameters
        # Boolean do_sex_check = false
        Boolean filter_discrepant_sex = true
        # Boolean do_sex_check_and_filter_discrepant_sex = (do_sex_check && filter_discrepant_sex)
        Float max_female_f = 0.2
        Float min_male_f = 0.8

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

    # Get multipliers
    Int variant_init_multiplier = floor((variant_init_count / 75000000) + 1)
    Int sample_init_multiplier = floor((sample_init_count / 750) + 1)
    call UTILS.wc as sample_ancestry_init_count {
        input:
            input_file = ancestry_samples,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    Int sample_ancestry_multiplier = floor((sample_ancestry_init_count.num_lines / 7500) + 1)
    
    scatter(i in range(length(chrs))) {

        String chr = "~{chrs[i]}"

        # Get variant multiplier
        call UTILS.wc as variant_count_chr {
            input:
                input_file = complete_bims[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
        Int variant_multiplier_chr = floor((variant_count_chr.num_lines / 7500000) + 1)

        Int combined_init_multiplier_chr = sample_init_multiplier * variant_multiplier_chr
        Int combined_ancestry_multiplier_chr = sample_ancestry_multiplier * variant_multiplier_chr

        # Partition data by ancestry and apply SNP call rate filter
        call PLINK.make_bed as all_snps_chr_subset_ancestry{
            input:
                bed_in = complete_beds[i],
                bim_in = complete_bims[i],
                fam_in = complete_fams[i],
                output_basename = "~{output_basename}_chr~{chr}_snp_miss",
                keep_samples = ancestry_samples,
                geno = max_variant_missing_call_rate,
                allow_no_sex = true,
                cpu = 1,
                mem_gb = 2 * combined_init_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Set het haploids to missing
        call PLINK.make_bed as all_snps_chr_het_hap_to_missing{
            input:
                bed_in = all_snps_chr_subset_ancestry.bed_out,
                bim_in = all_snps_chr_subset_ancestry.bim_out,
                fam_in = all_snps_chr_subset_ancestry.fam_out,
                set_hh_missing = true,
                allow_no_sex = true,
                output_basename = "~{output_basename}_chr~{chr}_snp_miss_het_hap_missing",
                cpu = 1,
                mem_gb = 4 * combined_ancestry_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Get list of SNPs that pass QC
        call UTILS.cut as all_snps_chr_post_snp_qc_snps{
            input:
                input_file = all_snps_chr_het_hap_to_missing.bim_out,
                args = "-f 2",
                output_filename = "~{output_basename}_chr~{chr}_post_snp_qc_snps.txt",
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get list of SNPs that pass QC
    call UTILS.cat as all_snps_post_snp_qc_snps{
        input:
            input_files = all_snps_chr_post_snp_qc_snps.output_file,
            output_filename = "~{output_basename}_post_snp_qc_snps.txt",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get combined multiplier
    Int combined_ancestry_multiplier = sample_ancestry_multiplier * variant_init_multiplier

    # Extract ancestry variants that pass SNP QC from ref SNP set
    call PLINK.make_bed as ref_snps_post_snp_qc{
        input:
            bed_in = filtered_bed,
            bim_in = filtered_bim,
            fam_in = filtered_fam,
            extract = all_snps_post_snp_qc_snps.output_file,
            keep_samples = ancestry_samples,
            allow_no_sex = true,
            output_basename = "~{output_basename}_post_snp_qc_ref_snps",
            cpu = 1 * combined_ancestry_multiplier,
            mem_gb = 4 * combined_ancestry_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get samples to keep based on call rate based on autosomes
    call PLINK.make_bed as ref_snps_get_low_called_samples{
        input:
            bed_in = ref_snps_post_snp_qc.bed_out,
            bim_in = ref_snps_post_snp_qc.bim_out,
            fam_in = ref_snps_post_snp_qc.fam_out,
            autosome = true,
            mind = max_sample_missing_call_rate,
            allow_no_sex = true,
            output_basename = "~{output_basename}_post_snp_qc_ref_snps_sample_miss",
            cpu = floor((0.5 * combined_ancestry_multiplier) + 1),
            mem_gb = 3 * combined_ancestry_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Keep samples that pass call rate filter
    call PLINK.make_bed as ref_snps_filter_low_called_samples{
        input:
            bed_in = ref_snps_post_snp_qc.bed_out,
            bim_in = ref_snps_post_snp_qc.bim_out,
            fam_in = ref_snps_post_snp_qc.fam_out,
            keep_samples = ref_snps_get_low_called_samples.fam_out,
            allow_no_sex = true,
            output_basename = "~{output_basename}_post_snp_qc_ref_snps_sample_miss",
            cpu = floor((0.5 * combined_ancestry_multiplier) + 1),
            mem_gb = 3 * combined_ancestry_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get list of samples with excess homozygosity
    call PLINK.get_excess_homo_samples{
        input:
            bed_in = ref_snps_filter_low_called_samples.bed_out,
            bim_in = ref_snps_filter_low_called_samples.bim_out,
            fam_in = ref_snps_filter_low_called_samples.fam_out,
            output_basename = "~{output_basename}_post_snp_qc_ref_snps_sample_miss_het",
            min_he = min_sample_he,
            max_he = max_sample_he,
            cpu = floor((0.5 * combined_ancestry_multiplier) + 1),
            mem_gb = 5 * combined_ancestry_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Filter them out (if they exist)
    if(size(get_excess_homo_samples.excess_homos) > 0){
        call PLINK.make_bed as ref_snps_remove_excess_homos{
            input:
                bed_in = ref_snps_filter_low_called_samples.bed_out,
                bim_in = ref_snps_filter_low_called_samples.bim_out,
                fam_in = ref_snps_filter_low_called_samples.fam_out,
                remove_samples = get_excess_homo_samples.excess_homos,
                allow_no_sex = true,
                output_basename = "~{output_basename}_post_snp_qc_ref_snps_sample_miss_het",
                cpu = floor((0.5 * combined_ancestry_multiplier) + 1),
                mem_gb = 3 * combined_ancestry_multiplier,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    File ref_snps_post_sample_qc_bed = select_first([ref_snps_remove_excess_homos.bed_out, ref_snps_filter_low_called_samples.bed_out])
    File ref_snps_post_sample_qc_bim = select_first([ref_snps_remove_excess_homos.bim_out, ref_snps_filter_low_called_samples.bim_out])
    File ref_snps_post_sample_qc_fam = select_first([ref_snps_remove_excess_homos.fam_out, ref_snps_filter_low_called_samples.fam_out])

    # Relatedness workflow to identify related samples
    call REL.relatedness_wf{
        input:
            bed_in = ref_snps_post_sample_qc_bed,
            bim_in = ref_snps_post_sample_qc_bim,
            fam_in = ref_snps_post_sample_qc_fam,
            output_basename = "~{output_basename}_relatedness",
            num_king_splits = num_king_splits,
            degree = degree,
            num_pcs_to_analyze = kinship_pcs_to_analyze,
            ancestral_pca_loading_cutoff = ancestral_pca_loading_cutoff,
            max_kinship_snps = max_kinship_snps,
            min_kinship_snps = min_kinship_snps,
            ancestral_pca_loading_step_size = ancestral_pca_loading_step_size,
            max_ancestral_pca_loading_cutoff = max_ancestral_pca_loading_cutoff,
            hwe_pvalue = hwe_filter_pvalue,
            max_variant_missing_call_rate = max_variant_missing_call_rate,
            ld_exclude_regions = ld_exclude_regions,
            ld_type = ld_type,
            window_size = ld_window_size,
            step_size = ld_step_size,
            r2_threshold = ld_r2_threshold,
            min_ld_maf = ld_maf_cutoff,
            king_cpu_per_split = 1 * combined_ancestry_multiplier,
            king_mem_gb_per_split = 4 * combined_ancestry_multiplier,
            pca_cpu = 1 * combined_ancestry_multiplier,
            pca_mem_gb = 16 * combined_ancestry_multiplier,
            qc_cpu = 1 * combined_ancestry_multiplier,
            qc_mem_gb = 4 * combined_ancestry_multiplier,
            ld_cpu = 1 * combined_ancestry_multiplier,
            ld_mem_gb = 4 * combined_ancestry_multiplier,
            merge_bed_cpu = 1 * combined_ancestry_multiplier,
            merge_bed_mem_gb = 4 * combined_ancestry_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Sex check to see if sex matches phenotype file
    call SEX.sex_check_wf as sex_check_wf{
        input:
            bed_in = ref_snps_post_sample_qc_bed,
            bim_in = ref_snps_post_sample_qc_bim,
            fam_in = ref_snps_post_sample_qc_fam,
            phenotype_file = ref_snps_post_sample_qc_fam,
            female_max_f = max_female_f,
            male_min_f = min_male_f,
            output_basename = "~{output_basename}_sex_check",
            header_rows = 0,
            fid_col = 0,
            iid_col = 1,
            sex_col = 4,
            delimiter = "tab",
            ld_exclude_regions = ld_exclude_regions,
            ld_type = ld_type,
            window_size = ld_window_size,
            step_size = ld_step_size,
            r2_threshold = ld_r2_threshold,
            min_ld_maf = ld_maf_cutoff,
            build_code = genome_build_code,
            no_fail = true,
            ld_cpu = 1 * combined_ancestry_multiplier,
            ld_mem_gb = 4 * combined_ancestry_multiplier,
            sex_check_cpu = 1 * combined_ancestry_multiplier,
            sex_check_mem_gb = 4 * combined_ancestry_multiplier,
            plink_cpu = 1 * combined_ancestry_multiplier,
            plink_mem_gb = 5 * combined_ancestry_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Optionally create filter lists from sex/relatedness workflows
    # Combine sample filter lists if filtering both
    # if(do_sex_check_and_filter_discrepant_sex && filter_related_samples){
    if(filter_discrepant_sex && filter_related_samples){
        call UTILS.get_file_union{
            input:
                input_files = [sex_check_wf.samples_to_remove, relatedness_wf.related_samples],
                output_filename = "~{output_basename}_sex_kin_remove",
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }
    # Set sex file if doing only sex removal
    # if(do_sex_check_and_filter_discrepant_sex && !filter_related_samples){
    if(filter_discrepant_sex && !filter_related_samples){
        File sex_remove = sex_check_wf.samples_to_remove
    }
    # Set kin file if doing only related removal
    # if(filter_related_samples && !do_sex_check_and_filter_discrepant_sex){
    if(filter_related_samples && !filter_discrepant_sex){
        File kin_remove = relatedness_wf.related_samples
    }
    File final_filter = select_first([get_file_union.output_file, sex_remove, kin_remove, write_lines([])])

    # Filter SNPs for HWE after identifying related samples to remove
    #### WHY ARE WE DOING THIS HERE? ####
    # Great question! For two reasons.
    # For 1, related individuals can mess up HWE calculations
    # For 2, we also want to use plink's '--nonfounders' option because otherwise it excludes any site with a parent in the fam file
    # Now, we can just run filter HWE sites based on all samples we KNOW are unrelated and completely ignore the parent/child/founder/nonfounder issues

    scatter(i in range(length(chrs))) {

        String hwe_chr = chrs[i]
        call HWE.hwe_filter_chr_wf as all_snps_chr_hwe_filter{
            input:
                bed_in = all_snps_chr_het_hap_to_missing.bed_out[i],
                bim_in = all_snps_chr_het_hap_to_missing.bim_out[i],
                fam_in = all_snps_chr_het_hap_to_missing.fam_out[i],
                output_basename = "~{output_basename}_chr~{hwe_chr}_final_basic_qc",
                chr = "~{hwe_chr}",
                build_code = genome_build_code,
                hwe_filter_pvalue = hwe_filter_pvalue,
                nonfounders = true,
                related_samples = relatedness_wf.related_samples,
                keep_samples = ref_snps_post_sample_qc_fam,
                plink_filter_cpu = 1 * combined_ancestry_multiplier_chr[i],
                plink_filter_mem_gb = 4 * combined_ancestry_multiplier_chr[i],
                plink_chr_cpu = 1 * combined_ancestry_multiplier_chr[i],
                plink_chr_mem_gb = 4 * combined_ancestry_multiplier_chr[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Apply those discrepant sex and related filters
        # if(do_sex_check_and_filter_discrepant_sex || filter_related_samples){
        if(filter_discrepant_sex || filter_related_samples){
            if(size(final_filter) > 0){
                # Filter out samples as determined by user options
                call PLINK.make_bed as all_snps_chr_filter_sex_and_kin{
                    input:
                        bed_in = all_snps_chr_hwe_filter.bed_out,
                        bim_in = all_snps_chr_hwe_filter.bim_out,
                        fam_in = all_snps_chr_hwe_filter.fam_out,
                        remove_samples = final_filter,
                        output_basename = "~{output_basename}_chr~{hwe_chr}_final_full_qc",
                        cpu = 1 * combined_ancestry_multiplier_chr[i],
                        mem_gb = 4 * combined_ancestry_multiplier_chr[i],
                        image_source = image_source,
                        ecr_repo = ecr_repo
                }
            }
        }

        # Pick final QC files
        File final_qc_merge_bed_by_chr = select_first([all_snps_chr_filter_sex_and_kin.bed_out, all_snps_chr_hwe_filter.bed_out])
        File final_qc_merge_bim_by_chr = select_first([all_snps_chr_filter_sex_and_kin.bim_out, all_snps_chr_hwe_filter.bim_out])
        File final_qc_merge_fam_by_chr = select_first([all_snps_chr_filter_sex_and_kin.fam_out, all_snps_chr_hwe_filter.fam_out])

    }
    
    # Merge chr bims
    Array[Pair[Array[File], String]] bim_sets_to_merge = [
        (complete_bims, "~{output_basename}_init.bim"),
        (all_snps_chr_subset_ancestry.bim_out, "~{output_basename}_post_snp_qc.bim"),
        (all_snps_chr_hwe_filter.bim_out, "~{output_basename}_final.bim")
    ]
    scatter(bims_to_merge in bim_sets_to_merge) {
        call UTILS.cat as merge_bims{
            input:
                input_files = bims_to_merge.left,
                output_filename = bims_to_merge.right,
                mem_gb = 2 * combined_ancestry_multiplier,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get SNP lists
    Array[Pair[File, String]] bims_for_snp_lists = [
        (merge_bims.output_file[0], "~{output_basename}_init_variants.txt"),
        (merge_bims.output_file[1], "~{output_basename}_passed_call_rate_variants.txt"),
        (merge_bims.output_file[2], "~{output_basename}_final_variants.txt")
    ]
    scatter(bim_for_snp_list in bims_for_snp_lists) {
        call UTILS.cut as snp_list{
            input:
                input_file = bim_for_snp_list.left,
                args =  "-f 2",
                output_filename = bim_for_snp_list.right,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get lists of removed SNPs
    Array[Pair[Array[File], String]] files_for_snp_comm = [
        ([snp_list.output_file[0], snp_list.output_file[1]], "~{output_basename}_low_call_rate_variants.txt"),
        ([snp_list.output_file[1], snp_list.output_file[0]], "~{output_basename}_failed_hwe_variants.txt")
    ]
    scatter(file_for_snp_comm in files_for_snp_comm) {
        call UTILS.comm as snp_comm{
            input:
                file1 = file_for_snp_comm.left[0],
                file2 = file_for_snp_comm.left[1],
                option = "-23",
                output_filename = file_for_snp_comm.right,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get SNP counts
    Array[File] files_for_snp_counts = [
        snp_list.output_file[0],
        snp_list.output_file[2],
        snp_comm.output_file[0],
        snp_comm.output_file[1]
    ]
    scatter(file_for_snp_count in files_for_snp_counts) {
        call UTILS.wc as snp_count {
            input:
                input_file = file_for_snp_count,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get lists of removed samples
    Array[Pair[Array[File], String]] files_for_sample_comm = [
        ([ref_snps_post_snp_qc.fam_out, ref_snps_filter_low_called_samples.fam_out], "~{output_basename}_samples_low_call_rate.fam"),
        ([ref_snps_filter_low_called_samples.fam_out, ref_snps_post_sample_qc_fam], "~{output_basename}_samples_excess_homozygosity.fam")
    ]
    scatter(file_for_sample_comm in files_for_sample_comm) {
        call UTILS.comm as sample_comm{
            input:
                file1 = file_for_sample_comm.left[0],
                file2 = file_for_sample_comm.left[1],
                option = "-23",
                output_filename = file_for_sample_comm.right,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get sample lists
    Array[Pair[File, String]] fams_for_sample_lists = [
        (ref_snps_post_snp_qc.fam_out, "~{output_basename}_samples_initial.tsv"),
        (final_qc_merge_fam_by_chr[0], "~{output_basename}_samples_final_full_qc.tsv"),
        (all_snps_chr_hwe_filter.fam_out[0], "~{output_basename}_samples_final_basic_qc.tsv"),
        (sample_comm.output_file[0], "~{output_basename}_samples_low_call_rate.tsv"),
        (sample_comm.output_file[1], "~{output_basename}_samples_excess_homozygosity.tsv")
    ]
    scatter(fam_for_sample_list in fams_for_sample_lists) {
        call UTILS.cut as sample_list{
            input:
                input_file = fam_for_sample_list.left,
                args =  "-f 1,2",
                output_filename = fam_for_sample_list.right,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Rename sample files
    Array[Pair[File, String]] sample_lists_to_rename = [
        (relatedness_wf.related_samples, "~{output_basename}_samples_related.tsv"),
        (sex_check_wf.samples_to_remove, "~{output_basename}_samples_failed_sex_check.tsv")
    ]
    scatter(sample_list_to_rename in sample_lists_to_rename) {
        call UTILS.rename_file as sample_rename{
            input:
                input_file = sample_list_to_rename.left,
                output_filename = sample_list_to_rename.right,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get samples counts
    Array[File] files_for_sample_counts = [
        sample_list.output_file[0],
        sample_list.output_file[1],
        sample_list.output_file[2],
        sample_list.output_file[3],
        sample_list.output_file[4],
        sample_rename.output_file[0],
        sample_rename.output_file[1]
    ]
    scatter(file_for_sample_count in files_for_sample_counts) {
        call UTILS.wc as sample_count {
            input:
                input_file = file_for_sample_count,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    output{

        # Fully filtered qc beds that have been through all filters
        FULL_QC_GENOTYPE_FILES full_qc_genotype_files = FULL_QC_GENOTYPE_FILES {
            beds: final_qc_merge_bed_by_chr,
            bims: final_qc_merge_bim_by_chr,
            fams: final_qc_merge_fam_by_chr
        }

        # QC files that have gone through basic site/sample filter but no kinship/sex filters
        # Can be used later if you don't want to use kinship/sex-discrepancy filtered datasets
        BASIC_QC_GENOTYPE_FILES basic_qc_genotype_files = BASIC_QC_GENOTYPE_FILES {
            beds: all_snps_chr_hwe_filter.bed_out,
            bims: all_snps_chr_hwe_filter.bim_out,
            fams: all_snps_chr_hwe_filter.fam_out
        }

        # Variant lists
        VARIANT_LISTS variant_lists = VARIANT_LISTS {
            initial: snp_list.output_file[0],
            final: snp_list.output_file[2],
            low_call_rate: snp_comm.output_file[0],
            failed_hwe: snp_comm.output_file[1]
        }

        # Sample lists
        SAMPLE_LISTS sample_lists = SAMPLE_LISTS {
            initial: sample_list.output_file[0],
            final_full_qc: sample_list.output_file[1],
            final_basic_qc: sample_list.output_file[2],
            low_call_rate: sample_list.output_file[3],
            excess_homozygosity: sample_list.output_file[4],
            related: sample_rename.output_file[0],
            failed_sex_check: sample_rename.output_file[1]
        }

        # Counts
        COUNTS counts = COUNTS {
            variant_initial: snp_count.num_lines[0],
            variant_final: snp_count.num_lines[1],
            variant_low_call_rate: snp_count.num_lines[2],
            variant_failed_hwe: snp_count.num_lines[3],
            sample_initial: sample_count.num_lines[0],
            sample_final_full_qc: sample_count.num_lines[1],
            sample_final_basic_qc: sample_count.num_lines[2],
            sample_low_call_rate: sample_count.num_lines[3],
            sample_excess_homozygosity: sample_count.num_lines[4],
            sample_related: sample_count.num_lines[5],
            sample_failed_sex_check: sample_count.num_lines[6]
        }

        # Relatedness wf outputs
        RELATEDNESS_WF_OUTPUTS relatedness_wf_outputs = RELATEDNESS_WF_OUTPUTS {
            annotated_kinships: relatedness_wf.annotated_kinship_output,
            kinship_id_maps: relatedness_wf.kinship_id_map
        }

        # Sex check wf outputs
        SEX_CHECK_WF_OUTPUTS sex_check_wf_outputs = SEX_CHECK_WF_OUTPUTS {
            sex_check_report: sex_check_wf.plink_sex_check_output
        }
    }
}
