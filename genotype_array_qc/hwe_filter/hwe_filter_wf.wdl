import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

workflow hwe_filter_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename
    Float hwe_filter_pvalue
    Int cpu
    Int mem_gb
    Int merge_cpu = 4
    Int merge_mem_gb

    # HWE filter on autosomes
    call PLINK.make_bed as hwe_auto{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            autosome = true,
            hwe_pvalue = hwe_filter_pvalue,
            output_basename = "${output_basename}.auto_only",
            cpu = cpu,
            mem_gb = mem_gb
    }

    # Merge X chr to ensure PAR/NONPAR are not split (chr 23 won't exist for split bed/bim/fam datasets)
    call PLINK.make_bed as merge_x_chr{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            output_basename = "${output_basename}.mergex",
            merge_x = true,
            merge_no_fail = true,
            cpu = cpu,
            mem_gb = mem_gb
    }

    # Check to see if dataset actually contains chrX
    call PLINK.contains_chr{
        input:
            bim_in = merge_x_chr.bim_out,
            chr = "23"
    }

    # Use females to compute HWE if sex chr present
    if(contains_chr.contains){

        # HWE filter on chrX for females
        call PLINK.make_bed as hwe_chrx{
            input:
                bed_in = merge_x_chr.bed_out,
                bim_in = merge_x_chr.bim_out,
                fam_in = merge_x_chr.fam_out,
                chr = 23,
                filter_females = true,
                hwe_pvalue = hwe_filter_pvalue,
                output_basename = "${output_basename}.chrX.females.hwe",
                cpu = cpu,
                mem_gb = mem_gb
        }

        # Get chrX SNPs in HWE
        call UTILS.cut{
            input:
                input_file = hwe_chrx.bim_out,
                args = "-f2",
                output_filename = "${output_basename}.chrX.females.extract"
        }

        # Extract chrX SNPs from full dataset
        call PLINK.make_bed as extract_hwe_chrx{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                extract = cut.output_file,
                output_basename = "${output_basename}.chrX.hwe",
                cpu = cpu,
                mem_gb = mem_gb
        }

        # Now autosome HWE snps with chrX hwe snps
        call PLINK.merge_two_beds as merge_auto_chrx{
            input:
                bed_in_a = hwe_auto.bed_out,
                bim_in_a = hwe_auto.bim_out,
                fam_in_a = hwe_auto.fam_out,
                bed_in_b = extract_hwe_chrx.bed_out,
                bim_in_b = extract_hwe_chrx.bim_out,
                fam_in_b = extract_hwe_chrx.fam_out,
                ignore_errors = false,
                output_basename = "${output_basename}",
                cpu = merge_cpu,
                mem_gb = merge_mem_gb
        }
    }

    output{
        File bed_out = select_first([merge_auto_chrx.bed_out, hwe_auto.bed_out])
        File bim_out = select_first([merge_auto_chrx.bim_out, hwe_auto.bim_out])
        File fam_out = select_first([merge_auto_chrx.fam_out, hwe_auto.fam_out])
    }

}