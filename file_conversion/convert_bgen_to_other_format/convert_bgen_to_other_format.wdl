import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow convert_bgen_to_other_format{
    Array[File] bgen_files
    Array[File] sample_files
    Array[File] output_basenames
    String export_format
    String ref_alt_mode
    Boolean rm_dup = false
    String rm_dup_mode = 'error'
    File? keep
    File? remove
    File? extract
    File? exclude

    Int cpu
    Int mem_gb

    # Parallelize
    scatter(i in range(length(bgen_files))){
        # Convert
        call PLINK.export_bgen_to_other_format{
            input:
                bgen_in = bgen_files[i],
                sample_in = sample_files[i],
                output_basename = output_basenames[i],
                export_format = export_format,
                ref_alt_mode = ref_alt_mode,
                rm_dup = rm_dup,
                rm_dup_mode = rm_dup_mode,
                keep = keep,
                remove = remove,
                extract = extract,
                exclude = exclude,
                cpu = cpu,
                mem_gb = mem_gb
        }

    }
}