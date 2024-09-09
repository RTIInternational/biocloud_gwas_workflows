version 1.1

import "plink.wdl" as PLINK
import "flashpca.wdl" as PCA
import "ld_prune_wf.wdl" as LD
import "king_kinship_wf.wdl" as KING

task restore_pedigree_ids{

    input {

        File id_list_in
        File id_map
        String output_filename

        # Runtime environment
        String docker_image = "rtibiocloud/tsv-utils:v2.2.0_5141a72"
        String ecr_image = "rtibiocloud/tsv-utils:v2.2.0_5141a72"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 1

    }

    command <<<
        set -e

        # Get number of records in file without pedigree
        input_id_lines=$(wc -l ~{id_list_in} | cut -d" " -f1)

        # Subset original ids in id_map based on ids in no_ped_lines
        tsv-join --filter-file ~{id_list_in} \
            --key-fields 1,2 \
            --delimiter " " \
            ~{id_map} | \
            cut -d" " -f3,4 > ~{output_filename}

        # Count number of records in mapped output file
        mapped_lines=$(wc -l ~{output_filename} | cut -d" " -f1)

        # Throw error if any ids didn't map
        if [ "$input_id_lines" -ne "$mapped_lines" ];
        then
            echo "Mapping error! Input id list contains $input_id_lines samples but mapped output contains $mapped_lines"
            exit 1
        fi
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File id_list_out = "~{output_filename}"
    }
}

workflow relatedness_wf{

    input {

        File bed_in
        File bim_in
        File fam_in
        String output_basename

        # QC Filtering cutoffs
        Float hwe_pvalue
        String? hwe_mode
        Float max_missing_site_rate
        Int qc_cpu = 1
        Int qc_mem_gb = 2

        # LD params
        File? ld_exclude_regions
        String ld_type
        Int window_size
        Int step_size
        Float r2_threshold
        Int ld_cpu = 1
        Int ld_mem_gb = 2
        Float min_ld_maf
        Int merge_bed_cpu = 1
        Int merge_bed_mem_gb = 4

        # King parameters
        Int degree = 3
        Int num_king_splits = 4
        Int king_cpu_per_split = 4
        Int king_mem_gb_per_split = 8

        # PCA parameters
        Int pca_cpu = 4
        Int pca_mem_gb = 8
        Int num_pcs_to_analyze = 3

        # Optionally specify genotype standardization method [binom | binom2]
        String? pca_standx
        # PCA Random seed
        Int? pca_seed
        # FlashPCA memory options (usually don't need to do anything unless very, very large data)
        Boolean? pca_batch
        Int? pca_blocksize

        # Ancestral SNP filtering parameters
        Float ancestral_pca_loading_cutoff = 0.003
        Int max_kinship_snps = 100000
        Int min_kinship_snps = 10000
        Float ancestral_pca_loading_step_size = 0.001
        Float max_ancestral_pca_loading_cutoff = 0.01

        # Runtime
        String image_source = "docker"
        String? ecr_repo

    }

    # Remove pedigree info from fam file
    call PLINK.remove_fam_pedigree{
        input:
            fam_in = fam_in,
            output_basename = "~{output_basename}.noped",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Do HWE/Call Rate/Filtering
    call PLINK.make_bed as init_qc_filter{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = remove_fam_pedigree.fam_out,
            output_basename = "~{output_basename}.qc",
            geno = max_missing_site_rate,
            hwe_pvalue = hwe_pvalue,
            hwe_mode = hwe_mode,
            cpu = qc_cpu,
            mem_gb = qc_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Do LD-prune of autosomes
    scatter(chr_index in range(22)){
        Int chr = chr_index + 1
        # Get subset of markers in LD
        call LD.ld_prune_wf as ld_prune{
            input:
                bed_in = init_qc_filter.bed_out,
                bim_in = init_qc_filter.bim_out,
                fam_in = init_qc_filter.fam_out,
                output_basename = "~{output_basename}.chr~{chr}.ldprune",
                ld_type = ld_type,
                window_size = window_size,
                step_size = step_size,
                r2_threshold = r2_threshold,
                cpu = ld_cpu,
                mem_gb = ld_mem_gb,
                maf = min_ld_maf,
                chr = chr,
                exclude_regions = ld_exclude_regions,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Merge chromosomes
    call PLINK.merge_beds{
        input:
            bed_in = ld_prune.bed_out,
            bim_in = ld_prune.bim_out,
            fam_in = ld_prune.fam_out,
            output_basename = "~{output_basename}.ldprune",
            cpu = merge_bed_cpu,
            mem_gb = merge_bed_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Call king to get related individuals to remove
    call KING.king_kinship_wf as round1_get_relateds{
        input:
            bed_in = merge_beds.bed_out,
            bim_in = merge_beds.bim_out,
            fam_in = merge_beds.fam_out,
            degree = degree,
            num_splits = num_king_splits,
            output_basename = "~{output_basename}.rd1.king",
            king_split_cpu = king_cpu_per_split,
            king_split_mem_gb = king_mem_gb_per_split,
            plink_cpu = qc_cpu,
            plink_mem_gb = qc_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo

    }

    if(size(round1_get_relateds.related_samples) > 0){

        # Remove related samples and redo QC
        call PLINK.make_bed as remove_round1_relateds{
            input:
                bed_in = merge_beds.bed_out,
                bim_in = merge_beds.bim_out,
                fam_in = merge_beds.fam_out,
                output_basename = "~{output_basename}.round1.unrelated",
                remove_samples = round1_get_relateds.related_samples,
                geno = max_missing_site_rate,
                hwe_pvalue = hwe_pvalue,
                hwe_mode = hwe_mode,
                cpu = qc_cpu,
                mem_gb = qc_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Re-do quality filter on unrelated samples
        call PLINK.make_bed as unrelated_qc_filter{
            input:
                bed_in = remove_round1_relateds.bed_out,
                bim_in = remove_round1_relateds.bim_out,
                fam_in = remove_round1_relateds.fam_out,
                output_basename = "~{output_basename}.round1.unrelated.qc",
                geno = max_missing_site_rate,
                hwe_pvalue = hwe_pvalue,
                hwe_mode = hwe_mode,
                cpu = qc_cpu,
                mem_gb = qc_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Do LD-prune of autosomes
        scatter(chr_index in range(22)){
            Int chr_unrelated = chr_index + 1
            # Get subset of markers in LD
            call LD.ld_prune_wf as ld_prune_unrelated{
                input:
                    bed_in = unrelated_qc_filter.bed_out,
                    bim_in = unrelated_qc_filter.bim_out,
                    fam_in = unrelated_qc_filter.fam_out,
                    output_basename = "~{output_basename}.chr~{chr_unrelated}.round1.unrelated.qc.ldprune",
                    ld_type = ld_type,
                    window_size = window_size,
                    step_size = step_size,
                    r2_threshold = r2_threshold,
                    cpu = ld_cpu,
                    mem_gb = ld_mem_gb,
                    maf = min_ld_maf,
                    chr = chr_unrelated,
                    exclude_regions = ld_exclude_regions,
                    image_source = image_source,
                    ecr_repo = ecr_repo
            }
        }

        # Merge chromosomes
        call PLINK.merge_beds as merge_ld_prune_unrelated{
            input:
                bed_in = ld_prune_unrelated.bed_out,
                bim_in = ld_prune_unrelated.bim_out,
                fam_in = ld_prune_unrelated.fam_out,
                output_basename = "~{output_basename}.round1.unrelated.qc.ldprune",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # PCA of remaining samples to identify ancestry-informative SNPs
    call PCA.flashpca as pca{
        input:
            bed_in = select_first([merge_ld_prune_unrelated.bed_out, merge_beds.bed_out]),
            bim_in = select_first([merge_ld_prune_unrelated.bim_out, merge_beds.bim_out]),
            fam_in = select_first([merge_ld_prune_unrelated.fam_out, merge_beds.fam_out]),
            ndim = num_pcs_to_analyze,
            standx = pca_standx,
            seed = pca_seed,
            cpu = pca_cpu,
            mem_gb = pca_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get list of non-ancestry informative SNPs from PC loadings
    call PCA.get_non_ancestry_informative_snps{
        input:
            pca_loadings = pca.loadings,
            output_filename = "~{output_basename}.nonacestry.snps.txt",
            loading_value_cutoff = ancestral_pca_loading_cutoff,
            max_snps = max_kinship_snps,
            min_snps = min_kinship_snps,
            cutoff_step_size = ancestral_pca_loading_step_size,
            max_cutoff = max_ancestral_pca_loading_cutoff,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Remove ancestry-informative SNPs from original QC dataset
    call PLINK.make_bed as remove_ancestry_snps{
        input:
            bed_in = merge_beds.bed_out,
            bim_in = merge_beds.bim_out,
            fam_in = merge_beds.fam_out,
            output_basename = "~{output_basename}.qc.nonancestral_snps",
            cpu = qc_cpu,
            mem_gb = qc_mem_gb,
            extract = get_non_ancestry_informative_snps.snps_to_keep,
            image_source = image_source,
            ecr_repo = ecr_repo
        }

    # Call king to get related individuals to remove
    call KING.king_kinship_wf as final_get_relateds{
        input:
            bed_in = remove_ancestry_snps.bed_out,
            bim_in = remove_ancestry_snps.bim_out,
            fam_in = remove_ancestry_snps.fam_out,
            degree = degree,
            num_splits = num_king_splits,
            output_basename = "~{output_basename}.final.king",
            king_split_cpu = king_cpu_per_split,
            king_split_mem_gb = king_mem_gb_per_split,
            plink_cpu = qc_cpu,
            plink_mem_gb = qc_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo

    }

    # Map related sample ids back to their original pedigree ids
    call restore_pedigree_ids{
        input:
            id_list_in = final_get_relateds.related_samples,
            id_map = remove_fam_pedigree.id_map_out,
            output_filename = "~{output_basename}.final.related_samples.remove",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        # Select either the merged kinship file (num_splits > 1) or the first shard (num_splits == 1)
        File related_samples = restore_pedigree_ids.id_list_out
        File annotated_kinship_output = final_get_relateds.annotated_kinship_output
        File nonancestral_snps = get_non_ancestry_informative_snps.snps_to_keep
        File kinship_id_map = remove_fam_pedigree.id_map_out
    }

}