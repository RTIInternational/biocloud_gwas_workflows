import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD
import "biocloud_gwas_workflows/biocloud_wdl_tools/king/king.wdl" as KING
import "biocloud_gwas_workflows/biocloud_wdl_tools/flashpca/flashpca.wdl" as PCA

workflow relatedness_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    # LD params
    File? ld_exclude_regions
    String ld_type
    Int window_size
    Int step_size
    Float r2_threshold
    Int ld_cpu = 1
    Int ld_mem_gb = 2
    Float min_ld_maf
    Int merge_bed_cpu = 4
    Int merge_bed_mem_gb = 8

    # Filtering cutoffs
    Float hwe_pvalue
    String? hwe_mode
    Float max_missing_site_rate
    Int qc_cpu = 4
    Int qc_mem_gb = 8

    # King parameters
    Int king_cpu = 16
    Int king_mem_gb = 32
    Int degree = 3

    # PCA parameters
    Int pca_cpu = 16
    Int pca_mem_gb = 16
    Int num_pcs_to_analyze = 3
    String? pca_standx
    Int? pca_div
    Int? pca_tol
    Boolean? pca_batch
    Int? pca_blocksize
    Int? pca_seed
    Int? pca_precision

    # Remove pedigree info from fam file
    call PLINK.remove_fam_pedigree{
        input:
            fam_in = fam_in,
            output_basename = "${output_basename}.noped"
    }

    # Do HWE/Call Rate/Filtering
    call PLINK.make_bed as init_qc_filter{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = remove_fam_pedigree.fam_out,
            output_basename = "${output_basename}.qc",
            geno = max_missing_site_rate,
            hwe_pvalue = hwe_pvalue,
            hwe_mode = hwe_mode,
            cpu = qc_cpu,
            cpu = qc_mem_gb
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
                output_basename = "${output_basename}.chr${chr}.ldprune",
                ld_type = ld_type,
                window_size = window_size,
                step_size = step_size,
                r2_threshold = r2_threshold,
                cpu = ld_cpu,
                mem_gb = ld_mem_gb,
                maf = min_ld_maf,
                chr = chr,
                exclude_regions = ld_exclude_regions
        }
    }

    # Merge chromosomes
    call PLINK.merge_beds{
        input:
            bed_in = ld_prune.bed_out,
            bim_in = ld_prune.bim_out,
            fam_in = ld_prune.fam_out,
            output_basename = "${output_basename}.ldprune",
            cpu = merge_bed_cpu,
            mem_gb = merge_bed_mem_gb
    }

    # Call king to get related individuals to remove
    call KING.king as king_unrelated{
        input:
            bed_in = merge_beds.bed_out,
            bim_in = merge_beds.bim_out,
            fam_in = merge_beds.fam_out,
            output_basename = "${output_basename}.king.",
            cpu = king_cpu,
            mem_gb = king_mem_gb,
            degree = degree,
            unrelated = true
    }

    if(size(king_unrelated.related_samples) > 0){
        # Convert related sample output file to plink-compatible sample list
        call KING.king_samples_to_ids as get_round1_relateds{
            input:
                king_samples_in = king_unrelated.related_samples,
                output_filename = "${output_basename}.rd1.related_samples.txt"
        }

        # Remove related samples and redo QC
        call PLINK.make_bed as remove_round1_relateds{
            input:
                bed_in = merge_beds.bed_out,
                bim_in = merge_beds.bim_out,
                fam_in = merge_beds.fam_out,
                output_basename = "${output_basename}.unrelated",
                cpu = qc_cpu,
                mem_gb = qc_mem_gb,
                remove_samples = get_round1_relateds.king_samples_out,
                geno = max_missing_site_rate,
                hwe_pvalue = hwe_pvalue,
                hwe_mode = hwe_mode,
                cpu = qc_cpu,
                cpu = qc_mem_gb
        }

        # Do LD-prune of autosomes
        scatter(chr_index in range(22)){
            Int chr_unrelated = chr_index + 1
            # Get subset of markers in LD
            call LD.ld_prune_wf as ld_prune_unrelated{
                input:
                    bed_in = remove_round1_relateds.bed_out,
                    bim_in = remove_round1_relateds.bim_out,
                    fam_in = remove_round1_relateds.fam_out,
                    output_basename = "${output_basename}.chr${chr_unrelated}.unrelated.ldprune",
                    ld_type = ld_type,
                    window_size = window_size,
                    step_size = step_size,
                    r2_threshold = r2_threshold,
                    cpu = ld_cpu,
                    mem_gb = ld_mem_gb,
                    maf = min_ld_maf,
                    chr = chr_unrelated,
                    exclude_regions = ld_exclude_regions
            }
        }

        # Merge chromosomes
        call PLINK.merge_beds as merge_ld_prune_unrelated{
            input:
                bed_in = ld_prune_unrelated.bed_out,
                bim_in = ld_prune_unrelated.bim_out,
                fam_in = ld_prune_unrelated.fam_out,
                output_basename = "${output_basename}.unrelated.ldprune",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb
        }
    }

    # PCA of remaining samples to identify any ancestral-level outliers
    call PCA.flashpca as flashpca{
        input:
            bed_in = select_first([merge_ld_prune_unrelated.bed_out, merge_beds.bed_out]),
            bim_in = select_first([merge_ld_prune_unrelated.bim_out, merge_beds.bim_out]),
            fam_in = select_first([merge_ld_prune_unrelated.fam_out, merge_beds.fam_out]),
            ndim = num_pcs_to_analyze,
            standx = pca_standx,
            div = pca_div,
            tol = pca_tol,
            seed = pca_seed,
            precision = pca_precision,
            cpu = pca_cpu,
            mem_gb = pca_mem_gb
    }

    output{
        File eigenvectors = flashpca.eigenvectors
        File pcs = flashpca.pcs
        File eigenvalues = flashpca.eigenvalues
        File pve = flashpca.pve
        File loadings = flashpca.loadings
        File meansd = flashpca.meansd
    }

}