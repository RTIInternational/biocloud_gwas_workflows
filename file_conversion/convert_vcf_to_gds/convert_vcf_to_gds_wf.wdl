import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_vcf_to_gds/convert_vcf_to_gds.wdl" as CONVERT

workflow convert_vcf_to_gds_wf{
    Array[File] input_files
    Array[File] output_files

    # Resources
    Int cpu = 1
    Int mem_gb = 4

    # Parallelize
    scatter(i in range(length(input_files))){
        # Convert files
        call CONVERT.convert_vcf_to_gds{
            input:
                in_file = input_files[i],
                out_file = output_files[i],
                cpu = cpu,
                mem_gb = mem_gb
        }

    }

    output{
        Array[File] results_files = output_files
    }
}