import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD
import "biocloud_gwas_workflows/biocloud_wdl_tools/terastructure/terastructure.wdl" as TSTRUCT
import "biocloud_gwas_workflows/biocloud_wdl_tools/terastructure_postprocessing/terastructure_postprocessing.wdl" as TSP

task get_structure_variants{
    File ref_bim
    File data_bim
    String output_filename
    Int num_snps = 10000

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 4
    Int mem_gb = 8

    command <<<
        set -e
        # Get lists of non-A/T and non-C/G SNPs
        perl -lane 'if (($F[4] eq "A" && $F[5] ne "T") || ($F[4] eq "T" && $F[5] ne "A") || ($F[4] eq "C" && $F[5] ne "G") || ($F[4] eq "G" && $F[5] ne "C")) { print $F[1]; }' \
            ${data_bim} | \
            sort -u | \
            grep "rs" \
            > data.variants

        # Get variants from full ref dataset
        cut -f 2,2 ${ref_bim} | \
            grep "rs" | \
            sort -u > ref.variants

        # Get intersection
        comm -12 data.variants ref.variants > intersect.variants

        # Get random subset
        perl -ne 'print rand()."\t".$_' intersect.variants | \
            sort -k1,1 | \
            head -${num_snps} | \
            cut -f2,2 > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File variant_ids = "${output_filename}"
    }

}

task get_ancestry_sample_ids{
    File ancestry_psam
    Array[String] ancestries_to_include
    Int sample_id_col
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2

    command <<<
        # Grep anything that matches any of the ancestries
        # And then select the sample IDs (adding 0 as the first column as a family id for plink1.9)
        grep "${sep="\\|" ancestries_to_include}" ${ancestry_psam} | awk '{print "0",$${sample_id_col}}' > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File sample_ids = "${output_filename}"
    }
}

workflow structure_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename
    Int max_snps_to_analyze = 10000

    # 1000G ref files
    File ref_bed
    File ref_bim
    File ref_fam

    # 1000G ancestry phenotype info
    String ancestry_psam
    Array[String] ancestries_to_include = ["EUR", "EAS", "AFR"]
    # Column where 1000G sample id is found in ancestry psam file
    Int ancestry_sample_id_col = 1
    # What level of ancestry is being examined [SUPERPOP | POP]
    String ancestry_pop_type = "SUPERPOP"
    # Cutoffs defining how to call ancestry from admixture proportions
    Array[String] ancestry_definitions = ["AFR=EAS<0.25;AFR>0.25", "EUR=EAS<0.25;AFR<0.25", "EAS=EAS>0.25;AFR<0.25"]

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

    # TeraStructur parameters
    Float terastructure_rfreq_perc = 0.2
    Int terastructure_cpu = 32
    Int terastructure_mem_gb = 32
    Boolean? force
    Int? seed
    Boolean? compute_beta
    File? id_map

    # Do LD-prune of autosomes
    scatter(chr_index in range(22)){
        Int chr = chr_index + 1
        # Get subset of markers in LD
        call LD.ld_prune_wf as ld_prune{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
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

    # Get list of SNPs to include
    call get_structure_variants{
        input:
            ref_bim = ref_bim,
            data_bim = merge_beds.bim_out,
            output_filename = "structure.variants",
            num_snps = max_snps_to_analyze
    }

    # Subset the selected overlapping SNPs from each dataset
    call PLINK.make_bed as subset_data{
        input:
            bed_in = merge_beds.bed_out,
            bim_in = merge_beds.bim_out,
            fam_in = merge_beds.fam_out,
            extract = get_structure_variants.variant_ids,
            snps_only = true,
            snps_only_type = "just-acgt",
            output_basename = "${output_basename}.structure_snps"
    }

    # Subset SNPs from ref dataset
    call get_ancestry_sample_ids{
        input:
            ancestry_psam = ancestry_psam,
            ancestries_to_include = ancestries_to_include,
            sample_id_col = ancestry_sample_id_col,
            output_filename = "${output_basename}.ancestry.samples"
    }

    # Get IDs of individuals of that ancestry
    call PLINK.make_bed as subset_ref{
        input:
            bed_in = ref_bed,
            bim_in = ref_bim,
            fam_in = ref_fam,
            keep_samples = get_ancestry_sample_ids.sample_ids,
            extract = get_structure_variants.variant_ids,
            snps_only = true,
            snps_only_type = 'just-acgt',
            output_basename = "${output_basename}.1000g.ref.structure_snps",
    }

    # Check to see if there are any merge conflicts that require strand-flipping
    call PLINK.merge_two_beds as get_merge_conflicts{
        input:
            bed_in_a = subset_data.bed_out,
            bim_in_a = subset_data.bim_out,
            fam_in_a = subset_data.fam_out,
            bed_in_b = subset_ref.bed_out,
            bim_in_b = subset_ref.bim_out,
            fam_in_b = subset_ref.fam_out,
            merge_mode = 7,
            ignore_errors = true,
            output_basename = "${output_basename}.combined.merge_conflicts"
    }

    # Flip ref SNPs if there are merge conflicts
    if(size(get_merge_conflicts.missnp_out) > 0){

        # Try flipping alleles for erroroneous SNPs
        call PLINK.make_bed as flip_ref{
            input:
                bed_in = subset_ref.bed_out,
                bim_in = subset_ref.bim_out,
                fam_in = subset_ref.fam_out,
                output_basename = "${output_basename}.1000g.ref.structure_snps.flipped",
                flip = get_merge_conflicts.missnp_out
        }
    }

    # Now actually try to do the merge
    call PLINK.merge_two_beds as combine_ref_and_data{
        input:
            bed_in_a = subset_data.bed_out,
            bim_in_a = subset_data.bim_out,
            fam_in_a = subset_data.fam_out,
            bed_in_b = select_first([flip_ref.bed_out, subset_ref.bed_out]),
            bim_in_b = select_first([flip_ref.bim_out, subset_ref.bim_out]),
            fam_in_b = select_first([flip_ref.fam_out, subset_ref.fam_out]),
            ignore_errors = false,
            output_basename = "${output_basename}.combined.structure_snps"
    }

    # Cluster dataset using terastructure
    call TSTRUCT.terastructure{
        input:
            bed_in = combine_ref_and_data.bed_out,
            bim_in = combine_ref_and_data.bim_out,
            fam_in = combine_ref_and_data.fam_out,
            k = length(ancestries_to_include),
            rfreq_perc = terastructure_rfreq_perc,
            force = force,
            seed = seed,
            compute_beta = compute_beta,
            id_map = id_map,
            cpu = terastructure_cpu,
            mem_gb = terastructure_mem_gb
    }

    # Analyze results and get ancestry groups from terastructure output
    call TSP.terastructure_postprocess{
        input:
            theta = terastructure.admixture_proportions,
            fam = combine_ref_and_data.fam_out,
            psam = ancestry_psam,
            ref_pop_type = ancestry_pop_type,
            ancestry_definitions = ancestry_definitions,
            output_basename = output_basename,
            dataset_label = "terrastructure_ancestry"
    }

    # Order keep files to be same order as input ancestry groups
    call TSP.order_by_ancestry{
        input:
            ancestry_files_in = terastructure_postprocess.ancestry_samples,
            ancestries = ancestries_to_include
    }

    output{
        Array[File] ancestry_thetas = terastructure_postprocess.ancestry_thetas
        Array[File] triangle_plots = terastructure_postprocess.triangle_plots
        Array[File] samples_by_ancestry = order_by_ancestry.ancestry_file_out
    }

}