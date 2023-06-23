import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD
import "biocloud_gwas_workflows/biocloud_wdl_tools/structure/structure.wdl" as STRUCT
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/structure_postprocessing/structure_postprocessing.wdl" as STRUCT_PP

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

        grep ${ancestry} $psam | awk '{print "0",$${sample_id_col}}' > ${output_filename}
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
    # Dataset containing samples you want to classify
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    ##### Define ancestry groups and 1000G reference samples you want to include in STRUCTURE run
    Boolean include_ref_samples = true
    File ref_bed
    File ref_bim
    File ref_fam
    # PSAM mappping ref_fam sample IDs to ancestry groups defined below (e.g. AFR, CEU)
    String ancestry_psam

    # Ancestries to pull ref samples from
    Array[String] ancestries_to_include = ["CEU", "CHB", "YRI"]
    # Ancestry pop type [SUPERPOP | POP] telling downstream modules which PSAM column to search for ancestry assignment
    String ancestry_pop_type = "POP"
    # Define how ancestry determinatoins are made from STRUCTURE admixture proportions
    # Must have 1 definition for each ancestry in ancestries_to_include
    Array[String] ancestry_definitions = ["YRI=CHB<0.25;YRI>0.25", "CEU=CHB<0.25;YRI<0.25", "CHB=CHB>0.25;YRI<0.25"]
    Array[String]? ancestry_aliases

    ##### STRUCTURE-specific options
    # Number of SNPs to subset
    Int max_snps_to_analyze = 10000
    # Whether or not to split samples into smaller groups for parallel processing (HIGHLY RECOMMENDED)
    Boolean split_mode = true
    # Max number of samples that will be included in a single STRUCTURE run when using split mode
    Int max_samples_per_split = 1000
    # Additional structure parameters that don't need to be changed in 99% of cases
    Int numreps = 1000
    Int burnin = 1000
    Int seed = 1523031945
    File? mainparams
    File? extraparams

    ##### LD-pruning parameters
    # Flag to determine whether structure SNPs should be LD-pruned
    Boolean do_ld_prune = true
    # Optinal list of high-linkage regions to exclude
    File? ld_exclude_regions
    String? ld_type
    Int? window_size
    Int? step_size
    Float? r2_threshold

    # Task-specific resources
    # STRUCTURE resources are per split so be careful if you're creating like 1000 splits. You don't actually need that much per job.
    Int structure_cpu = 8
    Int structure_mem_gb = 16
    Int ld_cpu = 1
    Int ld_mem_gb = 2
    Float? min_ld_maf
    Int merge_bed_cpu = 4
    Int merge_bed_mem_gb = 8
    Int plink_cpu = 1
    Int plink_mem_gb = 2
    Int structure_postprocess_cpu = 1
    Int structure_postprocess_mem_gb = 2

    # Set false to do unsuperviesed clustering (but realistically don't touch this)
    # Controls whether structure uses ref pop assignments, which it pretty much always should
    Boolean use_pop_info = true

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

    # Create STRUCTURE param files based on what output files will look like
    # Input will have no marker names and includes pop_data, pop_flag columns that will be used by STRUCTURE (USEPOPINFO=1)
    call STRUCT.make_structure_param_files{
        input:
            markernames = 0,
            pop_flag = 1,
            pop_data = 1,
            use_pop_info = if(use_pop_info) then 1 else 0,
            burnin = burnin,
            numreps = numreps,
            randomize = 0,
            pfrompopflagonly = 1
    }

    # Get list of SNPs to include
    call get_structure_variants{
        input:
            ref_bim = ref_bim,
            data_bim = sample_bim,
            output_filename = "structure.variants",
            num_snps = max_snps_to_analyze
    }

    # Get sample ids to extract from 1000G ref for each ancestry of interest
    scatter(ancestry_index in range(length(ancestries_to_include))){
        String ancestry = ancestries_to_include[ancestry_index]

        # Get list of samples from specified ancestry
        call get_ancestry_samples{
            input:
                ancestry_psam = ancestry_psam,
                ancestry = ancestry,
                sample_id_col = 1,
                output_filename = "${output_basename}.${ancestry}.samples"
        }

        # Assign a numerical pop id to these samples that will later be used by STRUCTURE
        call UTILS.append_column as add_pop_ids{
            input:
                input_file = get_ancestry_samples.sample_ids,
                value = ancestry_index + 1,
                output_filename = "${output_basename}.${ancestry}.samples.pop_id"
        }
    }

    # Cat sample ids into single list to be used to extract samples from reference dataset
    call UTILS.cat as cat_ancestry_ids{
        input:
            input_files = add_pop_ids.output_file,
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
            output_basename = basename(sample_fam, ".fam"),
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
                output_basename = "${output_basename}.split.${split_index}.data",
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
                output_basename = "${output_basename}.split.${split_index}.merge_conflicts",
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
                    output_basename = "${output_basename}.split.${split_index}.ref",
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
                output_basename = "${output_basename}.split.${split_index}.combined",
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
                output_basename = "${output_basename}.split.${split_index}.combined",
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }

        # Get actual number of snps to pass to structure
        call UTILS.wc as count_structure_snps{
            input:
                input_file = recode_to_ped.map_out
        }

        # Convert ped file to STRUCTURE dataset with pop information included for ref samples
        call STRUCT.ped2structure{
            input:
                ped_in = recode_to_ped.ped_out,
                ref_samples = cat_ancestry_ids.output_file,
                ref_pop_col = 2,
                ref_delim = "space",
                output_filename = "${output_basename}.split.${split_index}.structure.input"
        }

        # Cluster dataset using STRUCTURE
        call STRUCT.structure{
            input:
                mainparams=make_structure_param_files.mainparams_out,
                extraparams=make_structure_param_files.extraparams_out,
                input_file = ped2structure.structure_input,
                output_basename = "${output_basename}.split.${split_index}.structure",
                numloci = count_structure_snps.num_lines,
                k = length(ancestries_to_include),
                seed = seed,
                cpu = structure_cpu,
                mem_gb = structure_mem_gb
        }

        # Extract just admixture proportions and corresponding sample ids
        call STRUCT.parse_structure_output{
            input:
                structure_output = structure.structure_out,
                output_filename = "${output_basename}.split.${split_index}.structure.results.txt"
        }
    }

    # Combine FAM files into single file
    call UTILS.cat as merge_fam_files{
        input:
            input_files = combine_ref_and_data.fam_out,
            output_filename = "${output_basename}.structure.merged.fam"
    }

    # Combine structure ancestry proportions into single file
    call UTILS.cat as merge_structure_output{
        input:
            input_files = parse_structure_output.structure_out,
            output_filename = "${output_basename}.structure.results.merged.txt"
    }


    #Analyze results and get ancestry groups from terastructure output
    call STRUCT_PP.structure_postprocess{
        input:
            parsed_structure_output = merge_structure_output.output_file,
            fam = merge_fam_files.output_file,
            psam = ancestry_psam,
            ref_pop_type = ancestry_pop_type,
            ancestry_definitions = ancestry_definitions,
            output_basename = output_basename,
            cpu = structure_postprocess_cpu,
            mem_gb = structure_postprocess_mem_gb
    }

    # Order keep files to be same order as input ancestry groups
    call STRUCT_PP.order_by_ancestry{
        input:
            ancestry_files_in = structure_postprocess.ancestry_samples,
            ancestries = select_first([ancestry_aliases, ancestries_to_include])
    }

    # Count number of samples that weren't classified as anything by structure
    call UTILS.wc as count_unclassified_samples{
        input:
            input_file = structure_postprocess.unclassified_samples
    }

    output{
        Array[File] ancestry_proportions = structure_postprocess.ancestry_proportions
        Array[File] triangle_plots = structure_postprocess.triangle_plots
        Array[File] samples_by_ancestry = order_by_ancestry.ancestry_file_out
        Int unclassified_samples = count_unclassified_samples.num_lines
    }

}