import "biocloud_gwas_workflows/helper_workflows/cov_ldsc/cov_ldsc_reference_chr_wf.wdl" as COVLDSC_CHR

workflow cov_ldsc_wf{
    Array[File] all_beds
    Array[File] all_bims
    Array[File] all_fams

    File cov_file
    Array[String] chrs
    Array[String] output_file_prefixes

    
    scatter(index in range(length(chrs))){
        
        # Run cov-LDSC on that chromosome
        call COVLDSC_CHR.cov_ldsc_chr_wf as cov_ldsc_full{
            input:
                bed_in = all_beds[index],
                bim_in = all_bims[index],
                fam_in = all_fams[index],
                cov_file = cov_file,
                out_prefix = output_file_prefixes[index],
        }
    }


}
