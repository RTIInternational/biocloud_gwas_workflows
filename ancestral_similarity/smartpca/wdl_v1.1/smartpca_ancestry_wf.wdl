version 1.1

import "ld_prune_wf.wdl" as LD
import "assign_ancestry_mahalanobis.wdl" as ANCESTRY
import "convert_variant_ids.wdl" as IDCONVERT
import "plink.wdl" as PLINK
import "smartpca.wdl" as SMARTPCA
import "utils.wdl" as UTILS

workflow smartpca_ancestry_wf{

    input {

        ## Dataset containing samples you want to classify
        File dataset_bed
        File dataset_bim
        File dataset_fam
        String dataset_short_name
        String? dataset_display_name

        ## Reference 
        File ref_bed
        File ref_bim
        File ref_fam
        String ref_psam
        String ancestry_pop_type = "SUPERPOP"
        Array[String] ancestries_to_include = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        Array[String] ancestries_display_names = ["African", "American Admixed", "East Asian", "European", "South Asian"]
        String genome_build_code = "GRCh37"

        ## ID Conversion
        Boolean do_id_conversion = true
        Array[File] id_conversion_ref_files = []
        String? id_conversion_in_sep
        String? id_conversion_in_missing_allele
        String? id_conversion_in_deletion_allele
        String? id_conversion_ref_deletion_allele
        Int? id_conversion_in_chunk_size
        Int? id_conversion_ref_chunk_size
        Boolean? id_conversion_rescue_rsids

        ## LD pruning parameters
        String ld_type = "indep-pairwise"
        Int ld_window_size = 20000
        Int ld_step_size = 2000
        Float ld_maf_cutoff = 0.01
        Float ld_r2_threshold = 0.5
        File ld_exclude_regions

        ## Smartpca parameters
        String altnormstyle = "YES"
        Int numoutevec = 10
        Int numoutlieriter = 5

        ## Assign ancestry parameters
        Int pc_count = 10
        Int use_pcs_count = 10
        String midpoint_formula = "median"
        
        # Resources
        String image_source = "docker"
        String? ecr_repo

    }

    # Get sample counts for dataset and ref
    Array[File] files_for_sample_counts = [
        ref_fam,
        dataset_fam
    ]
    scatter(i in range(length(files_for_sample_counts))) {
        call UTILS.wc as sample_count {
            input:
                input_file = files_for_sample_counts[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Get variant counts for dataset and ref
    Array[File] files_for_snp_counts = [
        ref_bim,
        dataset_bim
    ]
    scatter(i in range(length(files_for_snp_counts))) {
        call UTILS.wc as snp_count {
            input:
                input_file = files_for_snp_counts[i],
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    Int ref_sample_multiplier = floor((sample_count.num_lines[0] / 2000) + 1)
    Int ref_snp_multiplier = floor((snp_count.num_lines[0] / 12000000) + 1)
    Int ref_mem_multiplier = ref_sample_multiplier * ref_snp_multiplier
    Int dataset_sample_multiplier = floor((sample_count.num_lines[1] / 2000) + 1)
    Int dataset_snp_multiplier = floor((snp_count.num_lines[1] / 12000000) + 1)
    Int dataset_mem_multiplier = dataset_sample_multiplier * dataset_snp_multiplier
    Int merged_mem_multiplier = ref_mem_multiplier + dataset_mem_multiplier

    # Get keep list of reference samples for specified ancestries
    call get_ref_samples{
        input:
            ancestry_psam = ref_psam,
            ancestry_pop_type = ancestry_pop_type,
            ancestries = ancestries_to_include,
            output_filename = "ref_samples.tsv",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Add pop ids to reference fam file
    call add_pop_ids_to_fam_files{
        input:
            dataset_fam_in = dataset_fam,
            ref_fam_in = ref_fam,
            ref_pop_xref = get_ref_samples.ref_samples,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    
    scatter(chr_index in range(22)){

        Int chr = chr_index + 1

        # Split ref by chromosome
        call PLINK.make_bed as split_ref_by_chr{
            input:
                bed_in = ref_bed,
                bim_in = ref_bim,
                fam_in= add_pop_ids_to_fam_files.ref_fam_out,
                keep_samples = get_ref_samples.ref_samples,
                snps_only = true,
                snps_only_type = 'just-acgt',
                chr = chr,
                output_basename = "ref_chr~{chr}",
                cpu = 2,
                mem_gb = 2 * ref_mem_multiplier,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Split dataset by chr
        call PLINK.make_bed as split_dataset_by_chr{
            input:
                bed_in = dataset_bed,
                bim_in = dataset_bim,
                fam_in= add_pop_ids_to_fam_files.dataset_fam_out,
                snps_only = true,
                snps_only_type = 'just-acgt',
                chr = chr,
                output_basename = "~{dataset_short_name}_chr~{chr}",
                cpu = 2,
                mem_gb = 2 * dataset_mem_multiplier,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        if(do_id_conversion){
            # Convert dataset variant IDs
            call IDCONVERT.convert_variant_ids as convert_dataset_variant_ids{
                input:
                    chr = chr,
                    in_file = split_dataset_by_chr.bim_out,
                    in_header = 0,
                    in_sep = "tab",
                    in_id_col = 1,
                    in_chr_col = 0,
                    in_pos_col = 3,
                    in_a1_col = 4,
                    in_a2_col = 5,
                    in_missing_allele = "0",
                    in_deletion_allele = "-",
                    in_chunk_size = 50000,
                    ref = id_conversion_ref_files[chr_index],
                    ref_deletion_allele = ".",
                    ref_chunk_size = 1000000,
                    output_filename = "~{dataset_short_name}_chr~{chr}_~{genome_build_code}_ids.bim",
                    cpu = 2,
                    mem_gb = 4,
                    image_source = image_source,
                    ecr_repo = ecr_repo
            }
        }

        File post_id_conversion_dataset_bim = select_first([convert_dataset_variant_ids.output_file, split_dataset_by_chr.bim_out])

        # Get overlapping variants between ref and datset
        call get_variants{
            input:
                ref_bim = split_ref_by_chr.bim_out,
                dataset_bim = post_id_conversion_dataset_bim,
                output_filename = "variants_chr~{chr}",
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Get variant counts for dataset and ref
        Array[File] files_for_snp_counts_chr = [
            split_ref_by_chr.bim_out,
            post_id_conversion_dataset_bim
        ]
        scatter(i in range(length(files_for_snp_counts))) {
            call UTILS.wc as snp_count_chr {
                input:
                    input_file = files_for_snp_counts_chr[i],
                    image_source = image_source,
                    ecr_repo = ecr_repo
            }
        }

        Int ref_sample_multiplier_chr = floor((sample_count.num_lines[0] / 2000) + 1)
        Int ref_snp_multiplier_chr = floor((snp_count_chr.num_lines[0] / 12000000) + 1)
        Int ref_mem_multiplier_chr = ref_sample_multiplier_chr * ref_snp_multiplier_chr
        Int dataset_sample_multiplier_chr = floor((sample_count.num_lines[1] / 2000) + 1)
        Int dataset_snp_multiplier_chr = floor((snp_count_chr.num_lines[1] / 12000000) + 1)
        Int dataset_mem_multiplier_chr = dataset_sample_multiplier_chr * dataset_snp_multiplier_chr
        Int merged_mem_multiplier_chr = ref_mem_multiplier_chr + dataset_mem_multiplier_chr

        # Subset SNPs and samples from ref dataset to get overlapping SNPs from desired ancestries
        call PLINK.make_bed as subset_ref{
            input:
                bed_in = split_ref_by_chr.bed_out,
                bim_in = split_ref_by_chr.bim_out,
                fam_in = split_ref_by_chr.fam_out,
                extract = get_variants.variant_ids,
                output_basename = "ref_chr~{chr}_smartpca_snps",
                cpu = 2,
                mem_gb = 1 * ref_mem_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Extract variants from dataset
        call PLINK.make_bed as subset_dataset{
            input:
                bed_in = split_dataset_by_chr.bed_out,
                bim_in = post_id_conversion_dataset_bim,
                fam_in = split_dataset_by_chr.fam_out,
                extract = get_variants.variant_ids,
                output_basename = "~{dataset_short_name}_chr~{chr}_smartpca_snps",
                cpu = 2,
                mem_gb = 1 * dataset_mem_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }


        # Check to see if there are any merge conflicts that require strand-flipping
        call PLINK.merge_two_beds as get_merge_conflicts{
            input:
                bed_in_a = subset_dataset.bed_out,
                bim_in_a = subset_dataset.bim_out,
                fam_in_a = subset_dataset.fam_out,
                bed_in_b = subset_ref.bed_out,
                bim_in_b = subset_ref.bim_out,
                fam_in_b = subset_ref.fam_out,
                merge_mode = 7,
                ignore_errors = true,
                output_basename = "ref_~{dataset_short_name}_chr~{chr}_merge_conflicts",
                cpu = 2,
                mem_gb = 1 * merged_mem_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Flip ref SNPs if there are merge conflicts
        if(size(get_merge_conflicts.missnp_out) > 0){

            # Try flipping alleles for erroroneous SNPs
            call PLINK.make_bed_plink1 as flip_ref{
                input:
                    bed_in = subset_ref.bed_out,
                    bim_in = subset_ref.bim_out,
                    fam_in = subset_ref.fam_out,
                    output_basename = "ref_flipped_chr~{chr}",
                    flip = get_merge_conflicts.missnp_out,
                    cpu = 2,
                    mem_gb = 1 * merged_mem_multiplier_chr,
                    image_source = image_source,
                    ecr_repo = ecr_repo
            }

        }

        # Now do the merge
        call PLINK.merge_two_beds as combine_ref_and_data{
            input:
                bed_in_a = subset_dataset.bed_out,
                bim_in_a = subset_dataset.bim_out,
                fam_in_a = subset_dataset.fam_out,
                bed_in_b = select_first([flip_ref.bed_out, subset_ref.bed_out]),
                bim_in_b = select_first([flip_ref.bim_out, subset_ref.bim_out]),
                fam_in_b = select_first([flip_ref.fam_out, subset_ref.fam_out]),
                ignore_errors = false,
                allow_no_sex = true,
                output_basename = "ref_~{dataset_short_name}_chr~{chr}",
                cpu = 2,
                mem_gb = 2 * merged_mem_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # LD prune
        call LD.ld_prune_wf as ld_prune{
            input:
                bed_in = combine_ref_and_data.bed_out,
                bim_in = combine_ref_and_data.bim_out,
                fam_in = combine_ref_and_data.fam_out,
                exclude_regions = ld_exclude_regions,
                output_basename = "ref_~{dataset_short_name}_chr~{chr}_ldpruned",
                ld_type = ld_type,
                window_size = ld_window_size,
                step_size = ld_step_size,
                r2_threshold = ld_r2_threshold,
                maf = ld_maf_cutoff,
                cpu = 4,
                mem_gb = 2 * merged_mem_multiplier_chr,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    # Merge LD pruned chromosomes into single dataset
    call PLINK.merge_beds as merge_chrs{
        input:
            bed_in = ld_prune.bed_out,
            bim_in = ld_prune.bim_out,
            fam_in = ld_prune.fam_out,
            allow_no_sex = true,
            output_basename = "ref_~{dataset_short_name}_ldpruned",
            cpu = 2,
            mem_gb = 1 * merged_mem_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Assign dummy IDs in bim and fam files to avoid smartpca error associated with long IDs; replace pop IDs with pop names
    call prepare_smartpca_input_files{
        input:
            bim_in = merge_chrs.bim_out,
            fam_in = merge_chrs.fam_out,
            pop_id_xref = get_ref_samples.ref_pop_id_xref,
            dataset_name = dataset_short_name,
            image_source = image_source,
            ecr_repo = ecr_repo
    }
    
    # Run smartpca
    call SMARTPCA.smartpca{
        input:
            genotypename = merge_chrs.bed_out,
            snpname = prepare_smartpca_input_files.bim_out,
            indivname = prepare_smartpca_input_files.fam_out,
            output_basename = dataset_short_name,
            altnormstyle = altnormstyle,
            numoutevec = numoutevec,
            numoutlieriter = numoutlieriter,
            poplist = ancestries_to_include,
            numthreads = 8,
            cpu = 8,
            mem_gb = 2 * merged_mem_multiplier,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    call process_smartpca_results{
        input:
            evec_in = smartpca.evec,
            eval_in = smartpca.eval,
            snpweight_in = smartpca.snpweight,
            fam_id_xref = prepare_smartpca_input_files.fam_id_xref,
            bim_id_xref = prepare_smartpca_input_files.bim_id_xref,
            output_basename = dataset_short_name,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    call ANCESTRY.assign_ancestry_mahalanobis {
        input:
            file_pcs = process_smartpca_results.evec_out,
            pc_count = pc_count,
            dataset = dataset_short_name,
            dataset_legend_label = dataset_display_name,
            ref_pops = ancestries_to_include,
            ref_pops_legend_labels = ancestries_display_names,
            use_pcs_count = use_pcs_count,
            midpoint_formula = midpoint_formula,
            cpu = 8,
            mem_gb = 8,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Order keep files to be same order as input ancestry groups
    call order_by_ancestry{
        input:
            ancestry_keep_files_in = assign_ancestry_mahalanobis.dataset_ancestry_keep_lists,
            ancestries = ancestries_to_include,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File evec = process_smartpca_results.evec_out
        File eval = process_smartpca_results.eval_out
        File snpweight = process_smartpca_results.snpweight_out
        File smartpca_log = smartpca.log
        Array[File] pre_processing_pc_plots  = assign_ancestry_mahalanobis.pre_processing_pc_plots
        File ref_dropped_samples  = assign_ancestry_mahalanobis.ref_dropped_samples
        File ref_raw_ancestry_assignments = assign_ancestry_mahalanobis.ref_raw_ancestry_assignments
        File ref_raw_ancestry_assignments_summary = assign_ancestry_mahalanobis.ref_raw_ancestry_assignments_summary
        File dataset_ancestry_assignments = assign_ancestry_mahalanobis.dataset_ancestry_assignments
        File dataset_ancestry_assignments_summary = assign_ancestry_mahalanobis.dataset_ancestry_assignments_summary
        Array[File] dataset_ancestry_assignments_plots = assign_ancestry_mahalanobis.dataset_ancestry_assignments_plots
        Array[File] dataset_ancestry_outliers_plots = assign_ancestry_mahalanobis.dataset_ancestry_outliers_plots
        Array[File] dataset_ancestry_keep_lists = order_by_ancestry.ancestry_keep_files_out
        Array[String] ancestry_labels = ancestries_to_include
    }

}

task get_ref_samples{

    input {

        File ancestry_psam
        String ancestry_pop_type
        Array[String] ancestries
        String output_filename

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 2

    }
    
    command <<<
        set -e

        psam=~{ancestry_psam}

         # Unzip psam if necessary
        if [[ ~{ancestry_psam} =~ \.gz$ ]]; then
            gunzip -c ~{ancestry_psam} > ref.psam
            psam=ref.psam
        fi

        ancestry_col=$(head -n 1 ~{ancestry_psam} | perl -ane 'for($i=0; $i<@F; $i++) { if (uc($F[$i]) eq "'~{ancestry_pop_type}'") { print $i; last; }}')
        perl -lane '
            use warnings;
            BEGIN {
                %ancestries = ();
                $nextAncestryId = 2;
                open(ANCESTRIES, "'~{write_lines(ancestries)}'");
                while (<ANCESTRIES>) {
                    chomp;
                    $ancestries{$_} = $nextAncestryId++;
                }
                close ANCESTRIES;
                open(XREF, ">ancestry_id_xref.tsv");
                print XREF "$_\t$ancestries{$_}" for (keys %ancestries);
                close XREF;
            }
            if (exists($ancestries{$F['$ancestry_col']})) {
                print join("\t","0",$F[0],$ancestries{$F['$ancestry_col']});
            }
        ' $psam > ~{output_filename}

    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File ref_samples = "~{output_filename}"
        File ref_pop_id_xref = "ancestry_id_xref.tsv"
    }

}

task add_pop_ids_to_fam_files {

    input {

        File dataset_fam_in
        File ref_fam_in
        File ref_pop_xref

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 2

    }

    command <<<

        set -e
        mkdir input_files

         # Unzip dataset fam if necessary
        if [[ ~{dataset_fam_in} =~ \.gz$ ]]; then
            gunzip -c ~{dataset_fam_in} > input_files/dataset.fam
        else
            ln -s ~{dataset_fam_in} input_files/dataset.fam
        fi

         # Unzip ref fam if necessary
        if [[ ~{ref_fam_in} =~ \.gz$ ]]; then
            gunzip -c ~{ref_fam_in} > input_files/ref.fam
        else
            ln -s ~{ref_fam_in} input_files/ref.fam
        fi

        # Dataset
        perl -lane 'use warnings; print join("\t",@F[0..4],"1");' input_files/dataset.fam > dataset.fam

        # Reference
        perl -lane '
            use warnings;
            BEGIN {
                %pop_xref = ();
                open(POP_XREF, "'~{ref_pop_xref}'");
                while (<POP_XREF>) {
                    chomp;
                    @cols = split("\t");
                    $pop_xref{$cols[0]."_".$cols[1]} = $cols[2];
                }
                close POP_XREF;
            }
            if (exists($pop_xref{$F[0]."_".$F[1]})) {
                print join("\t", @F[0..4], $pop_xref{$F[0]."_".$F[1]});
            } else {
                print join("\t",@F);
            }
        ' input_files/ref.fam > ref.fam
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File dataset_fam_out = "dataset.fam"
        File ref_fam_out = "ref.fam"
    }

}

task get_variants{

    input {

        File ref_bim
        File dataset_bim
        String output_filename

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 2

    }

    command <<<
        set -e

        dataset_bim=~{dataset_bim}
        ref_bim=~{ref_bim}

        # Unzip bims if necessary
        if [[ ~{dataset_bim} =~ \.gz$ ]]; then
            gunzip -c ~{dataset_bim} > data.bim
            dataset_bim=data.bim
        fi

        if [[ ~{ref_bim} =~ \.gz$ ]]; then
            gunzip -c ~{ref_bim} > ref.bim
            ref_bim=ref.bim
        fi

        # Get lists of non-A/T and non-C/G SNPs
        perl -lane 'use warnings; if (($F[4] eq "A" && $F[5] ne "T") || ($F[4] eq "T" && $F[5] ne "A") || ($F[4] eq "C" && $F[5] ne "G") || ($F[4] eq "G" && $F[5] ne "C")) { print $F[1]; }' \
            $dataset_bim | \
            sort -u | \
            grep "rs" \
            > data.variants

        # Get variants from full ref dataset
        cut -f 2,2 $ref_bim | \
            grep "rs" | \
            sort -u > ref.variants

        # Get intersection
        comm -12 data.variants ref.variants > ~{output_filename}

    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File variant_ids = "~{output_filename}"
    }

}

task prepare_smartpca_input_files {

    input {

        File bim_in
        File fam_in
        File pop_id_xref
        String dataset_name
        
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem = 2

    }

    command <<<

        # bim file
        perl -lane '
            use warnings;
            BEGIN {
                open(XREF, ">bim_id_xref.tsv");
                $nextId = 1;
            }
            $id = "ID_".$nextId++;
            print XREF join("\t",$id,$F[1]);
            print join("\t",$F[0],$id,@F[2..5]);
            END {
                close XREF;
            }
        ' ~{bim_in} > smartpca.bim

        # fam file
        perl -lane '
            use warnings;
            BEGIN {
                %popIdXref = ();
                $popIdXref{"1"} = "'~{dataset_name}'";
                open(POP_ID_XREF, "'~{pop_id_xref}'");
                while(<POP_ID_XREF>) {
                    chomp;
                    @cols = split;
                    $popIdXref{$cols[1]} = $cols[0];
                }
                close POP_ID_XREF;
                open(XREF, ">fam_id_xref.tsv");
                $nextId = 1;
            }
            $id = "ID_".$nextId++;
            print XREF join("\t",$id,$F[0]."___".$F[1]);
            print join("\t",$id,$id,@F[2..4],$popIdXref{$F[5]});
            END {
                close XREF;
            }
        ' ~{fam_in} > smartpca.fam
        
    >>>

    output {
        File bim_out = "smartpca.bim"
        File bim_id_xref = "bim_id_xref.tsv"
        File fam_out = "smartpca.fam"
        File fam_id_xref = "fam_id_xref.tsv"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem} GB"
    }

    parameter_meta {
        bim_in: "Plink formatted bim file."
        fam_in: "Plink formatted fam file."
        cpu: "Number of CPUs for the image."
        memory: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Rename the variant IDs in the bim file and the sample IDs in the fam file. smartpca will throw an error if the IDs are too long. So we changed them to a generic ID for the smartpca analysis, then incorporate them back in later."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task process_smartpca_results {

    input {

        File evec_in
        File eval_in
        File snpweight_in
        File bim_id_xref
        File fam_id_xref
        String output_basename

        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem = 2

    }

    command <<<

        # evec file
        echo -e "FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\tPOP" > ~{output_basename}_evec.tsv
        tail -n +2 ~{evec_in} |
            perl -lane '
                use warnings;
                BEGIN {
                    %famIdXref = ();
                    open(FAM_ID_XREF, "'~{fam_id_xref}'");
                    while(<FAM_ID_XREF>) {
                        /^(\S+)\t(\S+)___(\S+)/;
                        $famIdXref{$1.":".$1} = $2."\t".$3;
                    }
                    close FAM_ID_XREF;
                }
                print(join("\t", $famIdXref{$F[0]}, @F[1..11]));
            ' >> ~{output_basename}_evec.tsv
        
        # eval file
        perl -lne '/(\S+)/; print $1;' ~{eval_in} > ~{output_basename}_eval.tsv

        # snpweight file
        echo -e "VARIANT_ID\tCHR\tPOS\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10" > ~{output_basename}_snpweight.tsv
        tail -n +2 ~{snpweight_in} |
            perl -lane '
                use warnings;
                BEGIN {
                    %bimIdXref = ();
                    open(BIM_ID_XREF, "'~{bim_id_xref}'");
                    while(<BIM_ID_XREF>) {
                        /^(\S+)\t(\S+)/;
                        $bimIdXref{$1} = $2;
                    }
                    close BIM_ID_XREF;
                }
                print(join("\t", $bimIdXref{$F[0]}, @F[1..12]));
            ' >> ~{output_basename}_snpweight.tsv
    >>>

    output {
        File evec_out = "~{output_basename}_evec.tsv"
        File eval_out = "~{output_basename}_eval.tsv"
        File snpweight_out = "~{output_basename}_snpweight.tsv"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem} GB"
    }

}

task order_by_ancestry {

    input {

        # Utility for aligning ancestry keep files with a set of ancestry names
        Array[String] ancestry_keep_files_in
        Array[String] ancestries

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem = 2

    }

    command <<<
        set -e

        # Loop through each ancestry and output filenames with the ancestry prefix
        for ancestry in ~{sep=" " ancestries}; do
            grep -i $ancestry ~{write_lines(ancestry_keep_files_in)}
        done > "ordered_ancestry.txt"
    >>>

    output {
        Array[String] ancestry_keep_files_out = read_lines("ordered_ancestry.txt")
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem} GB"
    }

}
