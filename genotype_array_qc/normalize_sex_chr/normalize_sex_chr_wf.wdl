import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS


workflow normalize_sex_chr_wf{
    # Reformats a bed/bim/fam file if necessary to conform with a list of expected chromosomes
    # Conditionally splits/merges PAR/NONPAR regions of chrX based on list of expected chromosomes
    File bed_in
    File bim_in
    File fam_in

    # List of chromosomes you're expecting to be in bed/bim/fam
    Array[String] expected_chrs

    String output_basename
    String build_code
    Boolean no_fail = true

    Int plink_cpu
    Int plink_mem_gb

    # Check if we're expecting a chr23 in bim
    call UTILS.array_contains as expect_sex_chr{
        input:
            input_array = expected_chrs,
            query = "23"

    }

    # Check if we're expecting chr25 (NONPAR) in bim
    call UTILS.array_contains as expect_split_sex_chr{
        input:
            input_array = expected_chrs,
            query = "25"
    }

    # Check to see whether current bim sex chr is split into PAR/NONPAR
    call PLINK.contains_chr as bim_split_sex_chr{
        input:
            bim_in = bim_in,
            chr = "25"
    }

    # Split sex chr if it needs splitting
    # Case: You want PAR/NONPAR separated but bim only contains chr23
    if(expect_sex_chr.contains && expect_split_sex_chr.contains && !bim_split_sex_chr.contains){

        # Split X chr by PAR/NONPAR
        call PLINK.make_bed_plink1 as split_x_chr{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                output_basename = "${output_basename}.splitx",
                split_x = true,
                build_code = build_code,
                split_no_fail = no_fail,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }
    }

    # Merge sex chr if it needs merging
    # Case: You want PAR/NONPAR regions collapsed but bim contains 23 and 25
    if(expect_sex_chr.contains && !expect_split_sex_chr.contains && bim_split_sex_chr.contains){
        # Split X chr by PAR/NONPAR
        call PLINK.make_bed_plink1 as merge_x_chr{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                output_basename = "${output_basename}.mergex",
                merge_x = true,
                merge_no_fail = no_fail,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb
        }
    }

    output{
        # Return either the merge version,
        File bed_out = select_first([split_x_chr.bed_out, merge_x_chr.bed_out, bed_in])
        File bim_out = select_first([split_x_chr.bim_out, merge_x_chr.bim_out, bim_in])
        File fam_out = select_first([split_x_chr.fam_out, merge_x_chr.fam_out, fam_in])
    }
}