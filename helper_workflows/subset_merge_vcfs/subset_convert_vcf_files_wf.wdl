import "biocloud_gwas_workflows/helper_workflows/subset_merge_vcfs/subset_convert_vcf_files_chr_wf.wdl" as SUBSET_CONVERT

workflow subset_convert_vcf_wf{
    Array[File] vcfs_in
    File samples_file
    Array[String] output_filenames
    String output_type
    Array[String] chrs

    Array[String] output_basenames

    scatter(index in range(length(vcfs_in))){
    
        call SUBSET_CONVERT.convert_vcf_to_bed{
            input:
                vcf_in = vcfs_in[index],
                samples_file = samples_file,
                output_filename = output_filenames[index],
                output_type = output_type,
                output_basename = output_basenames[index]
        }
    }
    
    output {
        Array[File] bed = convert_vcf_to_bed.bed_out
        Array[File] bim = convert_vcf_to_bed.bim_out
        Array[File] fam = convert_vcf_to_bed.fam_out
        Array[File] log = convert_vcf_to_bed.plink_log
    }
    
}

