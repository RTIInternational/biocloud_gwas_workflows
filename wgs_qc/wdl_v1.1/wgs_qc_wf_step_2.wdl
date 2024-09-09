version 1.1

import "plink.wdl" as PLINK
import "hwe_filter_chr_wf.wdl" as HWE
import "relatedness_wf.wdl" as REL
import "sex_check_wf.wdl" as SEX
import "utils.wdl" as UTILS
# import "normalize_sex_chr_wf.wdl" as NORM

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

        # # Sample filter cutoffs
        Float max_sample_missing_rate = 0.03
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

        # General CPU/Mem for basically any PLINK task that runs on the whole dataset
        Int plink_cpu = 1
        Int plink_mem_gb = 2

        # General CPU/MEM for basically any PLINK task that runs on a single chr
        Int plink_chr_cpu = 1
        Int plink_chr_mem_gb = 2

        # General CPU/MEM for jobs that merge bed files either by chr or sample
        Int merge_bed_cpu = 4
        Int merge_bed_mem_gb = 8

        # Speicific tasks where resource limits may need to be adjusted for larger/smaller inputs
        Int sex_check_cpu = 4
        Int sex_check_mem_gb = 8
        Int king_cpu_per_split = 4
        Int king_mem_gb_per_split = 8
        Int pca_cpu = 4
        Int pca_mem_gb = 8

        # Container
        String image_source = "docker"
        String? ecr_repo
    
    }

    scatter(i in range(length(chrs))) {

        String chr = "~{chrs[i]}"

        # Partition data by ancestry and apply SNP call rate filter
        call PLINK.make_bed as all_snps_chr_subset_ancestry{
            input:
                bed_in = complete_beds[i],
                bim_in = complete_bims[i],
                fam_in = complete_fams[i],
                output_basename = "~{output_basename}_chr~{chr}_snp_miss",
                keep_samples = ancestry_samples,
                geno = max_missing_site_rate,
                allow_no_sex = true,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb,
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
                cpu = plink_cpu,
                mem_gb = plink_mem_gb,
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
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Filter samples based on call rate (autosomes)
    call PLINK.make_bed as ref_snps_filter_low_called_samples{
        input:
            bed_in = ref_snps_post_snp_qc.bed_out,
            bim_in = ref_snps_post_snp_qc.bim_out,
            fam_in = ref_snps_post_snp_qc.fam_out,
            autosome = true,
            mind = max_sample_missing_rate,
            allow_no_sex = true,
            output_basename = "~{output_basename}_post_snp_qc_ref_snps_sample_miss",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
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
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
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
                cpu = plink_cpu,
                mem_gb = plink_mem_gb,
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
            max_missing_site_rate = max_missing_site_rate,
            ld_exclude_regions = ld_exclude_regions,
            ld_type = ld_type,
            window_size = ld_window_size,
            step_size = ld_step_size,
            r2_threshold = ld_r2_threshold,
            min_ld_maf = ld_maf_cutoff,
            king_cpu_per_split = king_cpu_per_split,
            king_mem_gb_per_split = king_mem_gb_per_split,
            pca_cpu = pca_cpu,
            pca_mem_gb = pca_mem_gb,
            qc_cpu = plink_cpu,
            qc_mem_gb = plink_mem_gb,
            ld_cpu = plink_chr_cpu,
            ld_mem_gb = plink_chr_mem_gb,
            merge_bed_cpu = merge_bed_cpu,
            merge_bed_mem_gb = merge_bed_mem_gb,
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
            ld_cpu = plink_chr_cpu,
            ld_mem_gb = plink_chr_mem_gb,
            sex_check_cpu = sex_check_cpu,
            sex_check_mem_gb = sex_check_mem_gb,
            plink_cpu = plink_cpu,
            plink_mem_gb = plink_mem_gb,
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
                output_basename = "~{output_basename}_chr~{hwe_chr}_snp_miss_het_hap_missing_het_hwe",
                chr = "~{hwe_chr}",
                build_code = genome_build_code,
                hwe_filter_pvalue = hwe_filter_pvalue,
                nonfounders = true,
                related_samples = relatedness_wf.related_samples,
                keep_samples = ref_snps_post_sample_qc_fam,
                plink_filter_cpu = plink_cpu,
                plink_filter_mem_gb = plink_mem_gb,
                plink_chr_cpu = plink_chr_cpu,
                plink_chr_mem_gb = plink_chr_mem_gb,
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
                        output_basename = "~{output_basename}_chr~{hwe_chr}_post_snp_qc_post_sample_qc",
                        cpu = plink_cpu,
                        mem_gb = plink_mem_gb,
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
    
    # Merge bims for counts
    call UTILS.cat as merge_init_bims{
        input:
            input_files = complete_bims,
            output_filename = "~{output_basename}_init.bim",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.cat as merge_passed_snp_call_rate_bims{
        input:
            input_files = all_snps_chr_subset_ancestry.bim_out,
            output_filename = "~{output_basename}_post_snp_qc.bim",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.cat as merge_final_bims{
        input:
            input_files = all_snps_chr_hwe_filter.bim_out,
            output_filename = "~{output_basename}_final.bim",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get SNP lists
    call UTILS.cut as init_snp_list{
        input:
            input_file = merge_init_bims.output_file,
            args =  "-f 2",
            output_filename = "~{output_basename}_init_variants.txt",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.cut as passed_snp_call_rate_list{
        input:
            input_file = merge_passed_snp_call_rate_bims.output_file,
            args =  "-f 2",
            output_filename = "~{output_basename}_passed_call_rate_variants.txt",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.cut as final_snp_list{
        input:
            input_file = merge_final_bims.output_file,
            args =  "-f 2",
            output_filename = "~{output_basename}_final_variants.txt",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get lists of removed SNPs
    call UTILS.comm as low_call_rate_snps{
        input:
            file1 = init_snp_list.output_file,
            file2 = passed_snp_call_rate_list.output_file,
            option = "-23",
            output_filename = "~{output_basename}_low_call_rate_variants.txt",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.comm as failed_hwe_snps{
        input:
            file1 = passed_snp_call_rate_list.output_file,
            file2 = final_snp_list.output_file,
            option = "-23",
            output_filename = "~{output_basename}_failed_hwe_variants.txt",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get lists of removed samples
    call UTILS.comm as low_call_rate_samples{
        input:
            file1 = filtered_fam,
            file2 = ref_snps_filter_low_called_samples.fam_out,
            option = "-23",
            output_filename = "~{output_basename}_low_call_rate_samples.tsv",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    call UTILS.comm as excess_homozygosity_samples{
        input:
            file1 = ref_snps_filter_low_called_samples.fam_out,
            file2 = ref_snps_post_sample_qc_fam,
            option = "-23",
            output_filename = "~{output_basename}_excess_homozygosity_samples.tsv",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }


    output{

        # Fully filtered qc beds that have been through all filters
        Array[File] final_qc_beds_by_chr = final_qc_merge_bed_by_chr
        Array[File] final_qc_bims_by_chr = final_qc_merge_bim_by_chr
        Array[File] final_qc_fams_by_chr = final_qc_merge_fam_by_chr

        # QC files that have gone through basic site/sample filter but no kinship/sex filters
        # Can be used later if you don't want to use kinship/sex-discrepancy filtered datasets
        Array[File] basic_qc_beds_by_chr = all_snps_chr_hwe_filter.bed_out
        Array[File] basic_qc_bims_by_chr = all_snps_chr_hwe_filter.bim_out
        Array[File] basic_qc_fams_by_chr = all_snps_chr_hwe_filter.fam_out

        # Intermediate files for counts
        File variants_initial = init_snp_list.output_file
        File variants_final_bim = final_snp_list.output_file
        File variants_low_call_rate = low_call_rate_snps.output_file
        File variants_failed_hwe = failed_hwe_snps.output_file
        File samples_initial = filtered_fam
        File samples_final_basic = all_snps_chr_hwe_filter.fam_out[0]
        File samples_final = final_qc_merge_fam_by_chr[0]
        File samples_low_call_rate = low_call_rate_samples.output_file
        File samples_excess_homozygosity = excess_homozygosity_samples.output_file
        File samples_related = relatedness_wf.related_samples
        File samples_failed_sex_check = sex_check_wf.samples_to_remove

        # Kinship estimation outputs for each ancestry
        File annotated_kinships = relatedness_wf.annotated_kinship_output
        # Kinship id maps for each ancestry bc annotated kinship results have dummy family ids
        File kinship_id_maps =  relatedness_wf.kinship_id_map

        # Sex discrepancy outputs for each ancestry
        File sex_check_report= sex_check_wf.plink_sex_check_output
    }
}
