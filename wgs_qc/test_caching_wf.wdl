import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow test_caching_wf{
    File vcf
    String output_basename = "derp"
    Int cpu
    Int mem_gb
    String container_source = "docker"

    call PLINK.convert_vcf_to_bed as all_snps_chr_convert_vcf_to_bfile{
        input:
            vcf_in = vcf,
            output_basename = output_basename,
            cpu = cpu,
            mem_gb = mem_gb,
            container_source = container_source
    }

    output{
        File bed = all_snps_chr_convert_vcf_to_bfile.bed_out
        File bim = all_snps_chr_convert_vcf_to_bfile.bim_out
        File fam = all_snps_chr_convert_vcf_to_bfile.fam_out
    }
}

