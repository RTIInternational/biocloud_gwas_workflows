import "biocloud_gwas_workflows/helper_workflows/subset_merge_vcfs/subset_merge_vcf_files_chr_wf.wdl" as SUBSETMERGE

workflow subset_merge_wf{
    # Note this workflow is specific to input 5 studies which will be subsetted and merged
    Array[File] input_files
    Array[File] samples_files
    Array[String] output_filenames
    String maf_filter
    String output_type
    Array[String] merge_file_output_filenames
    
    # Resources
    Int cpu = 4
    Int mem_gb = 8

    # Parallelize
    scatter(i in range(length(input_files))){
        # Call per chr workflow
        call SUBSETMERGE.subset_merge_vcf_chr_wf{
            input:
                vcfs_in = input_files[i],
                samples_files = samples_files,
                output_filenames = output_filenames[i]
                maf_filter = maf_filter,
                output_type = output_type,
                merge_file_output_filename = merge_file_output_filenames[i]
        }

    }

}
