import "biocloud_gwas_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as SUBSET
import "biocloud_gwas_workflows/biocloud_wdl_tools/merge_utils/merge_utils.wdl" as MERGE

workflow subset_merge_vcf_chr_wf{
    # BCFtools view options
    Array[File] vcfs_in
    Array[File] samples_files
    Array[String] output_filenames
    String output_type
    String maf_filter
    ?String chr

    String merge_file_output_filename

    # Do subset workflow for each vcf file
    scatter(index in range(length(vcfs_in))){
        call SUBSET.view as view{
            input:
                vcf_in = vcfs_in[index],
                samples_file = samples_files[index],
                output_filename = output_filenames[index],
                maf_filter = maf_filter,
                output_type = output_type
        }
    }

    call MERGE.merge_vcfs as merge{
        input:
            input_vcfs = view.vcf_out,
            output_filename = merge_file_output_filename
    } 

    output {
        File merged_vcf = merge.merged_vcf
    }
}

