import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD
import "biocloud_gwas_workflows/biocloud_wdl_tools/terastructure/terastructure.wdl" as TSTRUCT
import "biocloud_gwas_workflows/biocloud_wdl_tools/terastructure_postprocessing/terastructure_postprocessing.wdl" as TSP
import "biocloud_gwas_workflows/biocloud_wdl_tools/terastructure_merger/terastructure_merger.wdl" as TSM
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

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

        data_bim=${data_bim}
        ref_bim=${ref_bim}

        # Unzip bims if necessary
        if [[ ${data_bim} =~ \.gz$ ]]; then
            gunzip -c ${data_bim} > data.bim
            data_bim=data.bim
        fi

        if [[ ${ref_bim} =~ \.gz$ ]]; then
            gunzip -c ${ref_bim} > ref.bim
            ref_bim=ref.bim
        fi

        # Get lists of non-A/T and non-C/G SNPs
        perl -lane 'if (($F[4] eq "A" && $F[5] ne "T") || ($F[4] eq "T" && $F[5] ne "A") || ($F[4] eq "C" && $F[5] ne "G") || ($F[4] eq "G" && $F[5] ne "C")) { print $F[1]; }' \
            $data_bim | \
            sort -u | \
            grep "rs" \
            > data.variants

        # Get variants from full ref dataset
        cut -f 2,2 $ref_bim | \
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
        set -e

        psam=${ancestry_psam}

         # Unzip psam if necessary
        if [[ ${ancestry_psam} =~ \.gz$ ]]; then
            gunzip -c ${ancestry_psam} > ref.psam
            psam=ref.psam
        fi

        grep "${sep="\\|" ancestries_to_include}" $psam | awk '{print "0",$${sample_id_col}}' > ${output_filename}
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

    Boolean split_mode = true
    Int max_samples_per_split = 2000

    # 1000G ref files
    Boolean include_ref_samples = true
    File ref_bed
    File ref_bim
    File ref_fam

    # 1000G ancestry phenotype info
    String ancestry_psam
    Array[String] ancestries_to_include = ["CEU", "CHB", "YRI"]
    # Column where 1000G sample id is found in ancestry psam file
    Int ancestry_sample_id_col = 1
    # What level of ancestry is being examined [SUPERPOP | POP]
    String ancestry_pop_type = "POP"
    # Cutoffs defining how to call ancestry from admixture proportions
    Array[String] ancestry_definitions = ["YRI=CHB<0.25;YRI>0.25", "CEU=CHB<0.25;YRI<0.25", "CHB=CHB>0.25;YRI<0.25"]

    # LD params
    Boolean do_ld_prune = true
    File? ld_exclude_regions
    String? ld_type
    Int? window_size
    Int? step_size
    Float? r2_threshold
    Int ld_cpu = 1
    Int ld_mem_gb = 2
    Float? min_ld_maf
    Int merge_bed_cpu = 4
    Int merge_bed_mem_gb = 8
    Int plink_cpu = 1
    Int plink_mem_gb = 2

    # TeraStructur parameters
    Float terastructure_rfreq_perc = 0.2
    Int terastructure_cpu = 32
    Int terastructure_mem_gb = 32
    Boolean? force
    Int? seed
    Boolean? compute_beta
    File? id_map

    # Optionally reduce to set of SNPs in LD
    if(do_ld_prune){
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

        # Merge LD chromosomes into single dataset
        call PLINK.merge_beds{
            input:
                bed_in = ld_prune.bed_out,
                bim_in = ld_prune.bim_out,
                fam_in = ld_prune.fam_out,
                output_basename = "${output_basename}.ldprune",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb
        }
    }

    File sample_bed = select_first([merge_beds.bed_out, bed_in])
    File sample_bim = select_first([merge_beds.bim_out, bim_in])
    File sample_fam = select_first([merge_beds.fam_out, fam_in])


    # Get list of SNPs to include
    call get_structure_variants{
        input:
            ref_bim = ref_bim,
            data_bim = sample_bim,
            output_filename = "structure.variants",
            num_snps = max_snps_to_analyze
    }

    # Get IDs of individuals from desired ancestry pops
    call get_ancestry_sample_ids{
        input:
            ancestry_psam = ancestry_psam,
            ancestries_to_include = ancestries_to_include,
            sample_id_col = ancestry_sample_id_col,
            output_filename = "${output_basename}.ancestry.samples"
    }

    # Subset SNPs and samples from ref dataset to get overlapping SNPs from desired ancestries
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
            cpu = plink_cpu,
            mem_gb = plink_mem_gb
    }

    # Split fam file into subsets of samples
    call SPLIT.split_file as get_sample_splits{
        input:
            input_file = sample_fam,
            output_basename = basename(sample_fam, ".fam")+".split",
            output_extension = ".fam",
            lines_per_split = max_samples_per_split,
    }

    # Do rest of structure WF in parallel for each batch of samples
    scatter(split_index in range(length(get_sample_splits.output_files))){
        # Subset the selected overlapping SNPs and samples from each split dataset
        call PLINK.make_bed as subset_data{
            input:
                bed_in = sample_bed,
                bim_in = sample_bim,
                fam_in = sample_fam,
                extract = get_structure_variants.variant_ids,
                keep_samples = get_sample_splits.output_files[split_index],
                snps_only = true,
                snps_only_type = "just-acgt",
                output_basename = "${output_basename}.structure_snps.split.${split_index}",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
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
                output_basename = "${output_basename}.combined.merge_conflicts.${split_index}",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb

        }

        # Flip ref SNPs if there are merge conflicts
        if(size(get_merge_conflicts.missnp_out) > 0){

            # Try flipping alleles for erroroneous SNPs
            call PLINK.make_bed as flip_ref{
                input:
                    bed_in = subset_ref.bed_out,
                    bim_in = subset_ref.bim_out,
                    fam_in = subset_ref.fam_out,
                    output_basename = "${output_basename}.1000g.ref.structure_snps.flipped.${split_index}",
                    flip = get_merge_conflicts.missnp_out,
                    cpu = plink_cpu,
                    mem_gb = plink_mem_gb
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
                output_basename = "${output_basename}.combined.structure_snps.${split_index}",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb
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

    }

    # Put terastructure theta split outputs into a standardized order
    # Otherwise the columns wouldn't necessarily align.
    File template_theta = terastructure.admixture_proportions[0]
    File template_fam   = combine_ref_and_data.fam_out[0]
    scatter(split_index in range(length(get_sample_splits.output_files))){
        call TSM.terastructure_merger{
            input:
                template_theta = template_theta,
                template_fam = template_fam,
                theta = terastructure.admixture_proportions[split_index],
                fam = combine_ref_and_data.fam_out[split_index],
                output_filename = "${output_basename}.${split_index}.standardized.theta.txt"
        }

    }

    # Combine FAM files into single file
    call UTILS.cat as merge_fam_files{
        input:
            input_files = combine_ref_and_data.fam_out,
            output_filename = "${output_basename}.combined.structure_snps.merged.fam"
    }

    # Combine theta files into single file
    call UTILS.cat as merge_theta_files{
        input:
            input_files = terastructure_merger.theta_out,
            output_filename = "${output_basename}.standardized.merged.theta.txt"
    }

    # Analyze results and get ancestry groups from terastructure output
    call TSP.terastructure_postprocess{
        input:
            theta = merge_theta_files.output_file,
            fam = merge_fam_files.output_file,
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

    # Count number of 1000G ref samples that weren't classified correctly
    call UTILS.wc as count_misclassified_ref_samples{
        input:
            input_file = terastructure_postprocess.misclassified_ref_samples
    }

    # Count number of samples that weren't classified as anything by structure
    call UTILS.wc as count_unclassified_samples{
        input:
            input_file = terastructure_postprocess.unclassified_samples
    }


    output{
        Array[File] ancestry_thetas = terastructure_postprocess.ancestry_thetas
        Array[File] triangle_plots = terastructure_postprocess.triangle_plots
        Array[File] samples_by_ancestry = order_by_ancestry.ancestry_file_out
        Int misclassified_ref_samples = count_misclassified_ref_samples.num_lines
        Int unclassified_samples = count_unclassified_samples.num_lines
    }

}