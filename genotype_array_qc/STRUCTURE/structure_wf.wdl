import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

task get_structure_variants{
    File ref_bim
    File data_bim
    String output_filename
    Int num_snps = 10000
    Int max_samples_per_structure_run = 10000

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
        comm -12 ${data_bim} ${ref_bim} > intersect.variants

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
    String ancestry
    Int sample_id_col
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 2

    command <<<
        grep "${ancestry}" ${ancestry_psam} | cut -f ${sample_id_col} > ${output_filename}
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
    Int num_snps_to_analyze = 10000

    # 1000G ref files
    File ref_bed
    File ref_bim
    File ref_fam

    # 1000G ancestry phenotype info
    String ancestry_psam
    Array[String] ancestries_to_include = ["EUR", "AMR", "EAS", "SAS", "AFR"]
    Int ancestry_sample_id_col = 1

    # LD params
    String ld_type
    Int window_size
    Int step_size
    Float r2_threshold
    Int ld_cpu
    Int ld_mem_gb
    Int min_ld_maf

    # Get Autosomes
    call PLINK.make_bed as get_autosome{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            autosome = true,
            output_basename = "${output_basename}.autosomes"
    }

    # Get subset of markers in LD
    call LD.ld_prune_wf as ld_prune{
        input:
            bed_in = get_autosome.bed_out,
            bim_in = get_autosome.bim_out,
            fam_in = get_autosome.fam_out,
            output_basename = "${output_basename}.ldprune",
            ld_type = ld_type,
            window_size = window_size,
            step_size = step_size,
            r2_threshold = r2_threshold,
            cpu = ld_cpu,
            mem_gb = ld_mem_gb,
            maf = min_ld_maf
    }

    # Get list of SNPs to include
    call get_structure_variants{
        input:
            ref_bim = ref_bim,
            data_bim = ld_prune.bim_out,
            output_filename = "structure.variants",
            num_snps = num_snps_to_analyze
    }

    # Subset the selected overlapping SNPs from each dataset
    call PLINK.make_bed as subset_data{
        input:
            bed_in = ld_prune.bed_out,
            bim_in = ld_prune.bim_out,
            fam_in = ld_prune.fam_out,
            extract = get_structure_variants.variant_ids,
            snps_only = true,
            snps_only_type = "just-acgt",
            output_basename = "${output_basename}.structure_snps"
    }

    # Creat STRUCTURE input file for each REF ancestry
    scatter(ancestry_index in range(length(ancestries_to_include))){
        String ancestry = ancestries_to_include[ancestry_index]
        # Get list of sample IDs for ancestry
        call get_ancestry_sample_ids{
            input:
                ancestry_psam = ancestry_psam,
                ancestry = ancestry,
                sample_id_col = ancestry_sample_id_col,
                output_filename = "{output_basename}.${ancestry}.samples"
        }

        # Get IDs of individuals of that ancestry
        call PLINK.make_bed as subset_ancestry{
            input:
                bed_in = ref_bed,
                bim_in = ref_bim,
                fam_in = ref_fam,
                keep = get_ancestry_sample_ids.sample_ids,
                extract = get_structure_variants.variant_ids,
                snps_only = true,
                snps_only_type = 'just-acgt',
                output_basename = "${output_basename}.1000g.${ancestry}.structure_snps",
        }

        # Convert to PED file
        call PLINK.bed_to_ped as ref_to_ped{
            input:
                bed_in = subset_ancestry.bed_out,
                bim_in = subset_ancestry.bim_out,
                fam_in = subset_ancestry.fam_out,
                output_basename = "${output_basename}.1000g.${ancestry}.structure_snps"
        }

        # Convert to STRUCTURE input format
        call PLINK.ped2structure as ref_to_structure{
            input:
                ped_in = ref_to_ped.ped_out,
                map_in = ref_to_ped.map_out,
                output_filename = "${output_basename}.1000g.${ancestry}.structure.input",
                pop_flag = 1,
                pop_to_use = ancestry_index
        }
    }

    # Cat Ref Structure files into single input file
    call UTILS.cat_files as cat_ref_structure{
        input:
            input_files = ref_to_structure.structure_input,
            output_filename = "${output_basename}.ref.structure.input"
    }

    # Split input into batches if large number of samples


    # Convert to PED format
    call PLINK.bed_to_ped as data_to_ped{
        input:
            bed_in = subset_data.bed_out,
            bim_in = subset_data.bim_out,
            fam_out = subset_data.fam_out,
            output_basename = "${output_basename}.structure_snsp"
    }



}