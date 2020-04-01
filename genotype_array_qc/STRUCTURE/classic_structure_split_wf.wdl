import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD
import "biocloud_gwas_workflows/biocloud_wdl_tools/structure/structure.wdl" as STRUCT
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT

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

task get_ancestry_samples{
    File ancestry_psam
    String ancestry
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

        grep "$ancestry" $psam | awk '{print "0",$${sample_id_col}}' > ${output_filename}
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
    Int max_samples_per_split = 1000

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
    Int structure_cpu = 32
    Int structure_mem_gb = 32
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

    # Get ids of 1000G samples separately for each ancestry
    scatter(ancestry in ancestries_to_include){
        call get_ancestry_samples{
            input:
                ancestry_psam = ancestry_psam,
                ancestry = ancestry,
                sample_id_col = ancestry_sample_id_col,
                output_filename = "${output_basename}.${ancestry}.samples"
        }
    }

    # Cat into single list of 1000G samples to extract from reference dataset
    call UTILS.cat as cat_ancestry_ids{
        input:
            input_files = get_ancestry_samples.sample_ids,
            output_filename = "${output_basename}.ancestry.samples"
    }

    # Subset SNPs and samples from ref dataset to get overlapping SNPs from desired ancestries
    call PLINK.make_bed as subset_ref{
        input:
            bed_in = ref_bed,
            bim_in = ref_bim,
            fam_in = ref_fam,
            keep_samples = cat_ancestry_ids.output_file,
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

        # Use split basename so all split filenames will appear in same order
        String split_basename = basename(get_sample_splits.output_files[split_index], ".fam")

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
                output_basename = "${split_basename}.data",
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
                output_basename = "${split_basename}.merge_conflicts",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb

        }

        # Flip ref SNPs if there are merge conflicts
        if(size(get_merge_conflicts.missnp_out) > 0){

            # Try flipping alleles for erroroneous SNPs
            call PLINK.make_bed_plink1 as flip_ref{
                input:
                    bed_in = subset_ref.bed_out,
                    bim_in = subset_ref.bim_out,
                    fam_in = subset_ref.fam_out,
                    output_basename = "${split_basename}.ref",
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
                output_basename = "${split_basename}.combined",
                cpu = merge_bed_cpu,
                mem_gb = merge_bed_mem_gb
        }

        ##### Convert to STRUCTURE input file
        # Convert to ped/map format
        call PLINK.recode_to_ped{
            input:
                bed_in = combine_ref_and_data.bed_out,
                bim_in = combine_ref_and_data.bim_out,
                fam_in = combine_ref_and_data.fam_out,
                output_basename = "${split_basename}.combined",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Get actual number of snps to pass to structure
        call UTILS.wc as count_structure_snps{
            input:
                input_file = recode_to_ped.map_out
        }

        # Convert to STRUCTURE dataset with pop information included for ref samples
        call STRUCT.ped2structure{
            input:
                ped_in = recode_to_ped.ped_out,
                pop_files = get_ancestry_samples.sample_ids,
                output_filename = "${split_basename}.combined.structure.input"
        }

        # Cluster dataset using terastructure
        call STRUCT.structure{
            input:
                mainparams = "",
                extraparams = "",
                stratparams = "",
                input_file = ped2structure.structure_input,
                output_basename = "${split_basename}.structure",
                numloci = count_structure_snps.num_lines,
                k = length(ancestries_to_include),
                seed = seed,
                cpu = structure_cpu,
                mem_gb = structure_mem_gb
        }

    }

    # Combine FAM files into single file
    call UTILS.cat as merge_fam_files{
        input:
            input_files = combine_ref_and_data.fam_out,
            output_filename = "${output_basename}.merged.fam"
    }

    # Need to compress into a tarzip for cat tsv files
    #call COLLECT.collect_large_file_list_wf as collect_thetas{
    #    input:
    #        input_files = structure.theta,
    #        output_dir_name = "${output_basename}_structure_output"
    #}

    # Combine theta files into single file
    #call TSV.append as merge_theta_files{
    #    input:
    #        input_files = collect_thetas.output_dir,
    #        output_filename = "${output_basename}.merged.theta.txt"
    #}

    # Analyze results and get ancestry groups from terastructure output
    #call TSP.terastructure_postprocess{
    #    input:
    #        theta = merge_theta_files.output_file,
    #        fam = merge_fam_files.output_file,
    #        psam = ancestry_psam,
    #        ref_pop_type = ancestry_pop_type,
    #        ancestry_definitions = ancestry_definitions,
    #        output_basename = output_basename,
    #        dataset_label = "terrastructure_ancestry"
    #}

    # Order keep files to be same order as input ancestry groups
    #call TSP.order_by_ancestry{
    #    input:
    #        ancestry_files_in = terastructure_postprocess.ancestry_samples,
    #        ancestries = ancestries_to_include
    #}

    # Count number of 1000G ref samples that weren't classified correctly
    #call UTILS.wc as count_misclassified_ref_samples{
    #    input:
    #        input_file = terastructure_postprocess.misclassified_ref_samples
    #}

    # Count number of samples that weren't classified as anything by structure
    #call UTILS.wc as count_unclassified_samples{
    #    input:
    #        input_file = terastructure_postprocess.unclassified_samples
    #}


    output{
        #Array[File] ancestry_thetas = terastructure_postprocess.ancestry_thetas
        #Array[File] triangle_plots = terastructure_postprocess.triangle_plots
        #Array[File] samples_by_ancestry = order_by_ancestry.ancestry_file_out
        #Int misclassified_ref_samples = count_misclassified_ref_samples.num_lines
        #Int unclassified_samples = count_unclassified_samples.num_lines
        Array[File] split_theta = structure.theta_out
    }

}