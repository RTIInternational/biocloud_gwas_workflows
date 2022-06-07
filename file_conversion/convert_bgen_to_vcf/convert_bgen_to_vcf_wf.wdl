import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow convert_bgen_to_vcf_wf{
    Array[File] bgen_files
    Array[File] sample_files
    Array[File] output_basenames
    String ref_alt_mode
    String vcf_dosage = "DS"
    File? keep
    File? remove
    File? extract
    File? exclude
    Boolean rm_dup = false
    String? rm_dup_mode

    Int cpu
    Int mem_gb

    # Parallelize
    scatter(i in range(length(bgen_files))){
        # Convert
        call PLINK.convert_bgen_to_vcf{
            input:
                bgen_in = bgen_files[i],
                sample_in = sample_files[i],
                output_basename = output_basenames[i],
                ref_alt_mode = ref_alt_mode,
                vcf_dosage = vcf_dosage,
                keep = keep,
                remove = remove,
                extract = extract,
                exclude = exclude,
                rm_dup = rm_dup,
                rm_dup_mode = rm_dup_mode,
                cpu = cpu,
                mem_gb = mem_gb
        }

    }

    output{
        Array[File] vcf_files = convert_bgen_to_vcf.vcf_out
        Array[File] log_files = convert_bgen_to_vcf.log_file
    }
}