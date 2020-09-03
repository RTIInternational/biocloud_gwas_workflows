import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_bgen_to_gds/convert_bgen_to_gds.wdl" as CONVERT

workflow convert_bgen_to_gds_wf{
    Array[File] input_files
    Array[File] output_files
    String storage_option = "LZMA_RA"
    String float_type = "double"
    Boolean geno = false
    Boolean dosage = false
    Boolean prob = false
    Boolean optimize = false
    Int parallel = 8

    # Resources
    Int mem_gb = 8

    # Parallelize
    scatter(i in range(length(input_files))){
        # Convert files
        call CONVERT.convert_bgen_to_gds{
            input:
                in_bgen = input_files[i],
                out_gds = output_files[i],
                storage_option = storage_option,
                float_type = float_type,
                geno = geno,
                dosage = dosage,
                prob = prob,
                optimize = optimize,
                parallel = parallel,
                mem_gb = mem_gb
        }

    }

    output{
        Array[File] results_files = output_files
    }
}