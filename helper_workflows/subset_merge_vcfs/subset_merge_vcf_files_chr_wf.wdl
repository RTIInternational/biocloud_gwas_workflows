import "biocloud_gwas_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as SUBSET
import "biocloud_gwas_workflows/biocloud_wdl_tools/bcftools_merge/bcftools_merge.wdl" as MERGE
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT

workflow subset_merge_vcf_chr_wf{
    # BCFtools view options
    Array[File] vcfs_in
    Array[File] samples_files
    Array[String] output_filenames
    String output_type
    Float maf_filter

    String file_out_prefix

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

    # Collect chunked sumstats files into single zip folder
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
            input_files = view.vcf_out,
            output_dir_name = file_out_prefix + "_extract_output"
    }

    output {
        File subset_vcfs = collect_sumstats.output_dir
    }
}

