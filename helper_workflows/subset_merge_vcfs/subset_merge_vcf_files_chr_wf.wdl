import "biocloud_gwas_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as SUBSET
import "biocloud_gwas_workflows/biocloud_wdl_tools/bcftools_merge/bcftools_merge.wdl" as MERGE

workflow subset_merge_vcf_chr_wf{
    # BCFtools view options
    Array[File] vcfs_in
    Array[File] samples_files
    Array[String] output_filenames
    String output_type
    Float maf_filter

    # BCFtools merge options
    File merge_files
    String output_filename_merged

    # Do subset workflow for each vcf file
    scatter(index in range(length(file_vcfs_in))){
        call SUBSET.view as view{
            input:
                vcf_in = vcfs_in[index],
                samples_file = samples_files[index],
                output_filename = output_filenames[index],
                maf_filter = maf_filter,
                output_type = output_type
        }
        #cat view.vcf_out >> merge_files
    }

    # Do merging of vcf files from the subset
    #call MERGE.merge as merge{
    #    input:
    #        file_in = merge_files,
    #        output_filename = output_filename_merged,
    #        output_type = output_type
    #}
    
 
    output {
        Array[File] subset_vcf = view.vcf_out
    }
}
