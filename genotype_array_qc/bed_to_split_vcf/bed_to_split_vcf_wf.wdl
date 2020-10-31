import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow bed_to_split_vcf_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    # Workflow can optionally be used on a single chr
    Array[String]? chrs

    Int plink_chr_cpu = 1
    Int plink_chr_mem_gb = 2

    # Get all chromosomes in bim if no list provided by user
    if(!defined(chrs)){
        call PLINK.get_bim_chrs{
            input:
                bim_in = bim_in
        }
    }

    Array[String] scatter_chrs = select_first([chrs, get_bim_chrs.chrs])

    scatter(chr in scatter_chrs){

        # Identify any samples that are missing most genotypes for a chr
        call PLINK.convert_bed_to_vcf as make_chr_vcf{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                chr = chr,
                output_basename = "${output_basename}.chr${chr}",
                cpu = plink_chr_cpu,
                mem_gb = plink_chr_mem_gb
        }
    }

    output{
        Array[File] vcf_out = make_chr_vcf.vcf_out
    }
}
