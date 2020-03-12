import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_to_impute2_ids/convert_to_impute2_ids.wdl" as IDCONVERT

workflow genotype_array_qc_wf{
    File bed
    File bim
    File fam
    Array[String] chrs
    Array[File] id_legend_files
    String output_basename

    # PAR/NONPAR Split/Merge parameters
    String build_code
    Boolean no_fail = true

    # Individual filtering options
    Float min_sample_call_rate

    # Resources for bed splitting/merging
    Int split_bed_cpu = 4
    Int split_bed_mem_gb = 8
    Int merge_bed_cpu = 2
    Int merge_bed_mem_gb = 12

    # Merge X chr to ensure PAR/NONPAR are not split (split-x will fail for pre-split files)
    call PLINK.make_bed as merge_x_chr{
        input:
            bed_in = bed,
            bim_in = bim,
            fam_in = fam,
            output_basename = "${output_basename}.mergex",
            merge_x = true,
            merge_no_fail = no_fail
    }

    # Split X chr by PAR/NONPAR
    call PLINK.make_bed as split_x_chr{
        input:
            bed_in = merge_x_chr.bed_out,
            bim_in = merge_x_chr.bim_out,
            fam_in = merge_x_chr.fam_out,
            output_basename = "${output_basename}.splitx",
            split_x = true,
            build_code = build_code,
            split_no_fail = no_fail
    }

    # Remove phenotype from fam file
    call PLINK.remove_fam_phenotype{
        input:
            fam_in = split_x_chr.fam_out,
            output_basename = "${output_basename}.splitx"
    }

    # Parallelize impute2 id conversion by chr
    scatter(chr_index in range(length(chrs))){
        String chr = chrs[chr_index]

        # Split by chr
        call PLINK.make_bed as split_bed{
            input:
                bed_in = split_x_chr.bed_out,
                bim_in = split_x_chr.bim_out,
                fam_in = remove_fam_phenotype.fam_out,
                output_basename = "${output_basename}.chr.${chr}",
                chr = chr,
                cpu = split_bed_cpu,
                mem_gb = split_bed_mem_gb
        }

        # Convert IDs to Impute2 format
        call IDCONVERT.convert_to_impute2_ids as impute2_id_bim{
            input:
                in_file = split_bed.bim_out,
                legend_file = id_legend_files[chr_index],
                contains_header = 0,
                id_col = 1,
                chr_col = 0,
                pos_col = 3,
                a1_col = 4,
                a2_col = 5,
                chr = chr,
                output_filename = basename(split_bed.bim_out),
                cpu = 2,
                mem_gb = 6
        }
    }

    # Merge into single file
    call PLINK.merge_beds{
        input:
            bed_in = split_bed.bed_out,
            bim_in = impute2_id_bim.output_file,
            fam_in = split_bed.fam_out,
            output_basename = "${output_basename}.impute2_id",
            cpu = merge_bed_cpu,
            mem_gb = merge_bed_mem_gb
    }

    # TODO: Remove duplicate IDS based on call rate

    # Remove failed subjects with >99% missing data
    call PLINK.make_bed as split_bed{
        input:
            bed_in = merge_beds.bed_out,
            bim_in = merge_beds.bim_out,
            fam_in = merge_beds.fam_out,
            output_basename = "${output_basename}.filter_failed_samples",
            mind = 0.99
    }

    # TODO: STRUCTURE WF to determine ancestry
    # TODO: Partition data by ancestry

    # For each ancestry group, do the following filtering steps
    # TODO: SNP Call rate filter
    # TODO: HWE filter
    # TODO: Set het haploids to missing
    # TODO: Subject call rate filter (autosomes)
    # TODO: Excess homozygosity filtering
    # TODO: Relatedness workflow
    # TODO: Sex check WF
    # TODO: Optionally Remove samples based on relatedness
    # TODO: Optionally Remove samples based on discrepant sex

    # TODO: Merge PAR/NONPAR regoins of chrX
    # TODO: Flag individuals missing chrX or other chr

    output{
        File bed_out = merge_beds.bed_out
        File bim_out = merge_beds.bim_out
        File fam_out = merge_beds.fam_out
        File fam_out_no_pheno = remove_fam_phenotype.fam_out
    }
}