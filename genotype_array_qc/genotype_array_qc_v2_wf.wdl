import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/id_conversion/convert_variant_ids_bfile/convert_variant_ids_bfile_wf.wdl" as IDCONVERT
import "biocloud_gwas_workflows/genotype_filtering/remove_duplicate_variants/remove_duplicate_variants_wf.wdl" as DUPLICATES
import "biocloud_gwas_workflows/ancestral_similarity/smartpca/smartpca_ancestry_wf.wdl" as ANCESTRY
import "biocloud_gwas_workflows/genotype_array_qc/hwe_filter/hwe_filter_wf.wdl" as HWE
import "biocloud_gwas_workflows/genotype_array_qc/relatedness/relatedness_wf.wdl" as REL
import "biocloud_gwas_workflows/genotype_array_qc/sex_check/sex_check_wf.wdl" as SEX
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

workflow genotype_array_qc_wf{

    # Input dataset parameters
    File bed
    File bim
    File fam
    String output_basename
    String dataset_short_name
    String dataset_display_name = dataset_short_name
    Array[String] chrs
    Array[String]? hwe_chrs # intended to be used in cases where the cohort is all-male or nearly all male and to tell the hwe_wf to skip chrx
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

    # Path to phenotype file containing sex calls
    # Intended to be the phenotype file that you will use for downstream association testing
    # If that's the fam file that's fine too but the important thing is to check that the sex calls that will
    # Be used in association testing are actually accurate
    File sex_check_phenotype_file
    # Number of header rows (all header rows must be removed)
    Int phenotype_header_rows
    # 0-based indices of required columns for sex check
    Int phenotype_fid_col
    Int phenotype_iid_col
    Int phenotype_sex_col
    String phenotype_delimiter

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

    # Remove phenotype from fam file
    call PLINK.remove_fam_phenotype{
        input:
            fam_in = fam,
            output_basename = "${output_basename}.no_pheno"
    }

    # Fix founder information in pedigree to make sure mother and father ids actually exist in dataset
    call PLINK.make_founders{
        input:
            fam_in = remove_fam_phenotype.fam_out,
            output_basename = "${output_basename}.no_pheno.make_founders"
    }

    # Count initial number of samples and sites
    call UTILS.wc as init_snp_count{
        input:
            input_file = bim
    }

    # Count initial number of samples and sites
    call UTILS.wc as init_sample_count{
        input:
            input_file = fam
    }

    # Convert variant IDs to standard format
    call IDCONVERT.convert_variant_ids_bfile_wf as convert_variant_ids{
        input:
            bed_in = bed,
            bim_in = bim,
            fam_in = make_founders.fam_out,
            ref_files = id_conversion_ref_files,
            chrs = chrs,
            output_basename = output_basename,
            in_sep = id_conversion_in_sep,
            in_missing_allele = id_conversion_in_missing_allele,
            in_deletion_allele = id_conversion_in_deletion_allele,
            ref_deletion_allele = id_conversion_ref_deletion_allele,
            in_chunk_size = id_conversion_in_chunk_size,
            ref_chunk_size = id_conversion_ref_chunk_size,
            rescue_rsids = id_conversion_rescue_rsids,
            convert_cpu = id_conversion_cpu,
            convert_mem_gb = id_conversion_mem_gb,
            plink_cpu = plink_cpu,
            plink_mem_gb = plink_mem_gb
    }

    # Remove duplicates
    call DUPLICATES.remove_duplicate_variants_wf as remove_duplicates{
        input:
            bed_in = convert_variant_ids.bed_out,
            bim_in = convert_variant_ids.bim_out,
            fam_in = convert_variant_ids.fam_out,
            output_basename = output_basename
    }

    # Count sites after removing dups
    call UTILS.wc as post_id_conversion_snp_count{
        input:
            input_file = remove_duplicates.bim_out
    }

    # Remove failed subjects with >99% missing data so they don't throw off anything downstream
    call PLINK.make_bed as filter_failed_samples{
        input:
            bed_in = remove_duplicates.bed_out,
            bim_in = remove_duplicates.bim_out,
            fam_in = remove_duplicates.fam_out,
            output_basename = "${output_basename}.filter_failed_samples",
            mind = 0.99,
            cpu = plink_cpu,
            mem_gb = plink_mem_gb
    }

    # Count samples after removing failed samples
    call UTILS.wc as failed_sample_count{
        input:
            input_file = filter_failed_samples.fam_out
    }

    # smartpca ancestry WF to partition by ancestry
    call ANCESTRY.smartpca_ancestry_wf{
        input:
            dataset_bed = filter_failed_samples.bed_out,
            dataset_bim = filter_failed_samples.bim_out,
            dataset_fam = filter_failed_samples.fam_out,
            dataset_short_name = dataset_short_name,
            dataset_display_name = dataset_display_name,
            ref_bed = ancestry_ref_bed,
            ref_bim = ancestry_ref_bim,
            ref_fam = ancestry_ref_fam,
            ref_psam = ancestry_ref_psam,
            ancestry_pop_type = ancestry_pop_type,
            ancestries_to_include = ancestries_to_include,
            ancestries_display_names = ancestries_display_names,
            ld_exclude_regions = ld_exclude_regions,
            plink_dataset_cpu = plink_cpu,
            plink_dataset_mem_gb = plink_mem_gb,
            convert_variant_ids_cpu = id_conversion_cpu,
            convert_variant_ids_mem_gb = id_conversion_mem_gb,
            smartpca_cpu = pca_cpu,
            smartpca_mem_gb = pca_mem_gb
    }

    # Filter out any ancestries with less than a minimum cutoff of samples
    scatter(ancestry_index in range(length(ancestries_to_include))){
        # Count number of samples in each ancestry
        call UTILS.wc as count_ancestry_samples{
            input:
                input_file = smartpca_ancestry_wf.dataset_ancestry_keep_lists[ancestry_index]
        }

        # Only include ancestries exceeding minimum number of samples
        if(count_ancestry_samples.num_lines > min_ancestry_samples_to_postprocess){
            Int ancestry_sample_count_maybe = count_ancestry_samples.num_lines
            String ancestries_maybe = ancestries_to_include[ancestry_index]
            File ancestry_samples_maybe = smartpca_ancestry_wf.dataset_ancestry_keep_lists[ancestry_index]
        }
    }

    # Remove null values from filtered ancestries to create list of ancestries to postprocess
    Array[Int] ancestry_sample_count = select_all(ancestry_sample_count_maybe)
    Array[String] postprocess_ancestries = select_all(ancestries_maybe)
    Array[File] postprocess_ancestry_samples = select_all(ancestry_samples_maybe)

    # Split by ancestry group and process each separately
    scatter(ancestry_index in range(length(postprocess_ancestries))){
        String ancestry = postprocess_ancestries[ancestry_index]
        File ancestry_samples = postprocess_ancestry_samples[ancestry_index]
        Int init_ancestry_sample_count = ancestry_sample_count[ancestry_index]

        # Partition data by ancestry and apply SNP call rate filter
        call PLINK.make_bed as subset_ancestry{
            input:
                bed_in = filter_failed_samples.bed_out,
                bim_in = filter_failed_samples.bim_out,
                fam_in = filter_failed_samples.fam_out,
                output_basename = "${output_basename}.${ancestry}.snp_miss",
                keep_samples = ancestry_samples,
                geno = max_missing_site_rate,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Count number of SNPs that didn't pass call rate filter
        call UTILS.wc as subset_ancestry_snp_count{
            input:
                input_file = subset_ancestry.bim_out
        }
        Int low_call_snp_count = post_id_conversion_snp_count.num_lines - subset_ancestry_snp_count.num_lines

        # Set het haploids to missing
        call PLINK.make_bed as het_hap_to_missing{
            input:
                bed_in = subset_ancestry.bed_out,
                bim_in = subset_ancestry.bim_out,
                fam_in = subset_ancestry.fam_out,
                set_hh_missing = true,
                output_basename = "${output_basename}.${ancestry}.snp_miss.het_hap_auto",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Get samples to filter based on call rate (autosomes)
        call PLINK.make_bed as get_low_called_samples{
            input:
                bed_in = het_hap_to_missing.bed_out,
                bim_in = het_hap_to_missing.bim_out,
                fam_in = het_hap_to_missing.fam_out,
                autosome = true,
                mind = max_sample_missing_rate,
                output_basename = "${output_basename}.${ancestry}.snp_miss.het_hap_miss",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Count samples removed due to low call rates
        call UTILS.wc as count_sample_call_filter{
            input:
                input_file = get_low_called_samples.fam_out
        }
        Int low_call_sample_count = init_ancestry_sample_count - count_sample_call_filter.num_lines

        # Filter out low call-rate samples
        call PLINK.make_bed as filter_low_called_samples{
            input:
                bed_in = het_hap_to_missing.bed_out,
                bim_in = het_hap_to_missing.bim_out,
                fam_in = het_hap_to_missing.fam_out,
                keep_samples = get_low_called_samples.fam_out,
                output_basename = "${output_basename}.${ancestry}.snp_miss.het_hap_miss.sample_miss",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Get list of samples with excess homozygosity
        call PLINK.get_excess_homo_samples{
            input:
                bed_in = filter_low_called_samples.bed_out,
                bim_in = filter_low_called_samples.bim_out,
                fam_in = filter_low_called_samples.fam_out,
                output_basename = "${output_basename}.${ancestry}.snp_miss.het_hap_miss.sample_miss.het",
                min_he = min_sample_he,
                max_he = max_sample_he,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Count number of samples removed due to excess homo
        call UTILS.wc as excess_homo_count{
            input:
                input_file = get_excess_homo_samples.excess_homos
        }
        Int excess_homo_sample_count = excess_homo_count.num_lines

        # Filter them out (if they exist)
        if(size(get_excess_homo_samples.excess_homos) > 0){
            call PLINK.make_bed as remove_excess_homos{
                input:
                    bed_in = filter_low_called_samples.bed_out,
                    bim_in = filter_low_called_samples.bim_out,
                    fam_in = filter_low_called_samples.fam_out,
                    remove_samples = get_excess_homo_samples.excess_homos,
                    output_basename = "${output_basename}.${ancestry}.snp_miss.hwe.het_hap_miss.sample_miss.het",
                    cpu = plink_cpu,
                    mem_gb = plink_mem_gb
            }
        }

        File qc_bed = select_first([remove_excess_homos.bed_out, filter_low_called_samples.bed_out])
        File qc_bim = select_first([remove_excess_homos.bim_out, filter_low_called_samples.bim_out])
        File qc_fam = select_first([remove_excess_homos.fam_out, filter_low_called_samples.fam_out])

        # Relatedness workflow to identify related samples
        call REL.relatedness_wf{
            input:
                bed_in = qc_bed,
                bim_in = qc_bim,
                fam_in = qc_fam,
                output_basename = "${output_basename}.${ancestry}",
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
        }

        # Count number of related samples that need to be removed
        call UTILS.wc as count_related_samples{
            input:
                input_file = relatedness_wf.related_samples
        }
        Int related_sample_removal_candidate_count = count_related_samples.num_lines

        # Filter SNPs for HWE after identifying related samples to remove
        #### WHY ARE WE DOING THIS HERE? ####
        # Great question! For two reasons.
        # For 1, related individuals can mess up HWE calculations
        # For 2, we also want to use plink's '--nonfounders' option because otherwise it excludes any site with a parent in the fam file
        # Now, we can just run filter HWE sites based on all samples we KNOW are unrelated and completely ignore the parent/child/founder/nonfounder issues
        call HWE.hwe_filter_wf{
            input:
                bed_in = qc_bed,
                bim_in = qc_bim,
                fam_in = qc_fam,
                chrs = hwe_chrs,
                hwe_filter_pvalue = hwe_filter_pvalue,
                nonfounders = true,
                related_samples = relatedness_wf.related_samples,
                output_basename = "${output_basename}.${ancestry}.snp_miss.het_hap_miss.sample_miss.het.hwe",
                plink_filter_cpu = plink_cpu,
                plink_filter_mem_gb = plink_mem_gb,
                plink_chr_cpu = plink_chr_cpu,
                plink_chr_mem_gb = plink_chr_mem_gb
        }

        # Count number of snps not in hwe
        call UTILS.wc as hwe_snp_count{
            input:
                input_file = hwe_filter_wf.bim_out
        }
        Int hwe_failed_snp_count = subset_ancestry_snp_count.num_lines - hwe_snp_count.num_lines

        # Sex check to see if sex matches phenotype file
        call SEX.sex_check_wf{
            input:
                bed_in = qc_bed,
                bim_in = qc_bim,
                fam_in = qc_fam,
                phenotype_file = sex_check_phenotype_file,
                female_max_f = max_female_f,
                male_min_f = min_male_f,
                output_basename = "${output_basename}.${ancestry}",
                header_rows = phenotype_header_rows,
                fid_col = phenotype_fid_col,
                iid_col = phenotype_iid_col,
                sex_col = phenotype_sex_col,
                delimiter = phenotype_delimiter,
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
                plink_mem_gb = plink_mem_gb
        }

        # Count number of sex-discrepant samples
        call UTILS.wc as count_sex_check_failed_samples{
            input:
                input_file = sex_check_wf.samples_to_remove
        }
        Int sex_check_failed_samples = count_sex_check_failed_samples.num_lines

        File hwe_qc_bed = hwe_filter_wf.bed_out
        File hwe_qc_bim = hwe_filter_wf.bim_out
        File hwe_qc_fam = hwe_filter_wf.fam_out

        # Optionally apply filter lists from sex/relatedness workflows to QC filtered samples
        # This section looks kinda weird but it's only because WDL doesn't have 'else/elseif'
        # Combine sample filter lists if filtering both
        if(filter_discrepant_sex && filter_related_samples){
            call UTILS.get_file_union{
                input:
                    input_files = [sex_check_wf.samples_to_remove, relatedness_wf.related_samples],
                    output_filename = "${output_basename}.sex.kin.remove"
            }
            String combined_suffix = "sex_check.unrelated"

            # Count length of samples in union
            call UTILS.wc as get_sex_kin_filter_count{
                input:
                    input_file = get_file_union.output_file
            }

            # Number of related samples that also failed sex check
            Int sex_related_overlap = sex_check_failed_samples + related_sample_removal_candidate_count - get_sex_kin_filter_count.num_lines
        }

        # Set sex file if doing only sex removal
        if(filter_discrepant_sex && !filter_related_samples){
            File sex_remove = sex_check_wf.samples_to_remove
            String sex_suffix = "sex_check"
        }

        # Set kin file if doing only related removal
        if(filter_related_samples && !filter_discrepant_sex){
            File kin_remove = relatedness_wf.related_samples
            String kin_suffix = "unrelated"
        }


        # Now actually apply those filters
        if(filter_discrepant_sex || filter_related_samples){
            # Select the file that contains the set of samples that need to be filtered
            # One of these will necessarily be defined
            File final_filter = select_first([get_file_union.output_file, sex_remove, kin_remove])
            File final_suffix = select_first([combined_suffix, sex_suffix, kin_suffix])

            # Remove samples if resulting file isn't empty
            if(size(final_filter) > 0){
                # Filter out samples as determined by user options
                call PLINK.make_bed as filter_sex_and_kin{
                    input:
                        bed_in = hwe_qc_bed,
                        bim_in = hwe_qc_bim,
                        fam_in = hwe_qc_fam,
                        remove_samples = final_filter,
                        output_basename = "${output_basename}.${ancestry}.snp_miss.hwe.het_hap_miss.sample_miss.het.${final_suffix}",
                        cpu = plink_cpu,
                        mem_gb = plink_mem_gb
                }
            }
        }

        # Get number of overlapping samples in related/sex check to account for differences
        # In cases of overlap, final_samples + related_samples + sex_failed_samples + homo_samples + low_call_samples will != init_ancestry samples
        # Need a convenient way to tell if someting is wrong so we correct explicitly for this overlap
        Int sex_related_overlap_count = select_first([sex_related_overlap, 0])


        # Finally re-merge PAR/NONPAR regions of chrX
        File split_final_bed = select_first([filter_sex_and_kin.bed_out, hwe_qc_bed])
        File split_final_bim = select_first([filter_sex_and_kin.bim_out, hwe_qc_bim])
        File split_final_fam = select_first([filter_sex_and_kin.fam_out, hwe_qc_fam])
        String final_basename = basename(split_final_bed, ".bed") + ".xmerged"

        call PLINK.make_bed as final_merge_x_chr{
            input:
                bed_in = split_final_bed,
                bim_in = split_final_bim,
                fam_in = split_final_fam,
                output_basename = final_basename,
                merge_x = true,
                merge_no_fail = true,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Final snp count
        call UTILS.wc as final_snp_count{
            input:
                input_file = split_final_bim
        }

        # Final sample count
        call UTILS.wc as final_sample_count{
            input:
                input_file = split_final_fam
        }

        # Check that number of samples removed makes sense
        Int total_removed_snps = init_snp_count.num_lines - final_snp_count.num_lines
        Int total_removed_samples = init_ancestry_sample_count - final_sample_count.num_lines
        Float pct_pass_snps = (final_snp_count.num_lines * 1.0)/(init_snp_count.num_lines * 1.0)
        Float pct_pass_samples = (final_sample_count.num_lines * 1.0)/(init_ancestry_sample_count * 1.0)

        # Check to make sure actual number of removed samples == expected removed samples
        # If this number is different from total_removed_samples something unexpected happened
        Int expected_removed_samples = sex_check_failed_samples + related_sample_removal_candidate_count + excess_homo_sample_count + low_call_sample_count - sex_related_overlap_count
        Int unaccounted_samples = total_removed_samples - expected_removed_samples
    }



    output{
        # Filter count metrics
        Int init_snps = init_snp_count.num_lines
        Int init_samples = init_sample_count.num_lines
        Int duplicate_snps = init_snp_count.num_lines - post_id_conversion_snp_count.num_lines
        Int failed_samples = init_sample_count.num_lines - failed_sample_count.num_lines
        Array[Int] samples_by_ancestry = count_ancestry_samples.num_lines
        Array[Int] low_call_snps_by_ancestry = low_call_snp_count
        Array[Int] hwe_failed_snps_by_ancestry = hwe_failed_snp_count
        Array[Int] low_call_samples_by_ancestry = low_call_sample_count

        Array[Int] excess_homo_samples_by_ancestry = excess_homo_sample_count
        Array[Int] related_samples_by_ancestry = related_sample_removal_candidate_count
        Array[Int] failed_sex_check_by_ancestry = sex_check_failed_samples
        Array[Int] final_qc_snps_by_ancestry = final_snp_count.num_lines
        Array[Int] final_qc_samples_by_ancestry = final_sample_count.num_lines
        Array[Int] total_removed_snps_by_ancestry = total_removed_snps
        Array[Int] total_removed_samples_by_ancestry = total_removed_samples
        Array[Float] pct_pass_snps_by_ancestry = pct_pass_snps
        Array[Float] pct_pass_samples_by_ancestry = pct_pass_samples
        Array[Int] unaccounted_samples_by_ancestry = unaccounted_samples

        # Fully filtered qc beds that have been through all filters
        Array[File] final_qc_bed = final_merge_x_chr.bed_out
        Array[File] final_qc_bim = final_merge_x_chr.bim_out
        Array[File] final_qc_fam = final_merge_x_chr.fam_out

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

        # Samples removed for low call rates or excess homozygosity for each ancestry
        Array[File] excess_homo_samples = get_excess_homo_samples.excess_homos

        # QC files that have gone through basic site/sample filter but no kinship/sex filters
        # Can be used later if you don't want to use kinship/sex-discrepancy filtered datasets
        Array[File] basic_qc_bed = hwe_qc_bed
        Array[File] basic_qc_bim = hwe_qc_bim
        Array[File] basic_qc_fam = hwe_qc_fam

        # Kinship estimation outputs for each ancestry
        Array[File] related_samples = relatedness_wf.related_samples
        Array[File] annotated_kinship = relatedness_wf.annotated_kinship_output
        # Kinship id maps for each ancestry bc annotated kinship results have dummy family ids
        Array[File] kinship_id_map =  relatedness_wf.kinship_id_map

        # Sex discrepancy outputs for each ancestry
        Array[File] sex_check_report = sex_check_wf.plink_sex_check_output
        Array[File] sex_check_failed_samples = sex_check_wf.samples_to_remove
    }
}
