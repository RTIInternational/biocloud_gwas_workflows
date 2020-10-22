import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_seq_gds_to_snp_gds/convert_seq_gds_to_snp_gds.wdl" as CONVERT

workflow convert_seq_gds_to_snp_gds_wf{
    Array[File] input_files
    Array[File] output_files
    String compress_geno = "ZIP_RA"
    Boolean optimize = false
    Boolean dosage = false

    # Resources
    Int mem_gb = 8
    Int cpu = 1

    # Parallelize
    scatter(i in range(length(input_files))){
        # Convert files
        call CONVERT.convert_seq_gds_to_snp_gds{
            input:
                in_seq_gds = input_files[i],
                out_snp_gds = output_files[i],
                compress_geno = compress_geno,
                optimize = optimize,
                dosage = dosage,
                cpu = cpu,
                mem_gb = mem_gb
        }

    }

    output{
        Array[File] results_files = output_files
    }
}