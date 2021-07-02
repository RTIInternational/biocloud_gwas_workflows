import "biocloud_gwas_workflows/helper_workflows/cov_ldsc/cov_ldsc_reference_chr_wf.wdl" as COVLDSC_CHR
import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plinl/wdl" as PLINK

workflow cov_ldsc_wf{
    File plink_format_file
    File = cov_file
    Array[String] = chrs
    String = output_file_prefix
    Int ldsc_cpu
    Int ldsc_mem_gb

    
    scatter(index in range(length(chrs))){
        
        # Split plink file by chromosome
        call PLINK.make_bed as make_bed{
            input:
                bed_in = plink_format_file.bed,
                bim_in = plink_format_file.bim,
                fam_in = plink_format_file.fam,
                output_basename = plink_format_file.byChr,
                chr = index,
        }

        # Run cov-LDSC on that chromosome
        call COVLDSC_CHR.cov_ldsc_chr_wf as cov_ldsc_full{
            input:
                plink_format_prefix = make_bed.bed_out,
                cov_file = cov_file,
                out_prefix = output_file_prefix.chr${index},
                ldsc_cpu = ldsc_cpu,
                ldsc_mem_gb = ldsb_mem_gb,
        }
    }


}
