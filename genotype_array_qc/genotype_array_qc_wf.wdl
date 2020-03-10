import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_to_impute2_ids/convert_to_impute2_ids.wdl" as IDCONVERT

workflow genotype_array_qc_wf{
    File bed
    File bim
    File fam
    Array[String] chrs
    Array[File] id_legend_files
    String output_basename

    # Resources for bed splitting/merging
    Int split_bed_cpu = 4
    Int split_bed_mem_gb = 8
    Int merge_bed_cpu = 2
    Int merge_bed_mem_gb = 12

    # Parallelize impute2 id conversion by chr
    scatter(chr_index in range(length(chrs))){
        String chr = chrs[chr_index]

        # Split by chr
        call PLINK.make_bed as split_bed{
            input:
                bed_in = bed,
                bim_in = bim,
                fam_in = fam,
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

    # Remove duplicate IDS based on call rate
    # Flag individuals missing chrX or other chr

    # Remove phenotype from fam file
    call PLINK.remove_fam_phenotype{
        input:
            fam_in = merge_beds.fam_out,
            output_basename = "${output_basename}.impute2_id.no_phenotype"
    }

    # Format phenotype file to standard format

    # Run STRUCTURE wf to determine ancestry

    # Partition data by ancestry

    # Call rate filter
    # HWE filter
    # Subject call rate filter (autosomes)
    # Relatedness workflow
    # Remove samples based on relatedness
    # Sex check and sample removal
    # Excess homozygosity filtering
    # Set het haploids to missing

    output{
        File bed_out = merge_beds.bed_out
        File bim_out = merge_beds.bim_out
        File fam_out = merge_beds.fam_out
        File fam_out_no_pheno = remove_fam_phenotype.fam_out
    }
}