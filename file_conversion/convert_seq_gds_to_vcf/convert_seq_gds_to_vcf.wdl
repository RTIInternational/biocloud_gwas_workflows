import "biocloud_gwas_workflows/biocloud_wdl_tools/seqarray_tools/seqarray_tools.wdl" as SEQARRAY

workflow convert_seq_gds_to_vcf{

    Array[File] input_gds_files
    Array[String] output_vcf_files
    Array[File] input_gds_variant_annot
    Array[File] input_gds_sample_annot
    File? variant_ids
    String? variant_annot_xref_col
    String? variant_annot_gds_id_col
    File? sample_ids
    String? sample_annot_xref_col
    String? sample_annot_gds_id_col
    Boolean? compress = true

    # Runtime attributes
    Int cpu = 1
    Int mem_gb = 4

    # Parallelize
    scatter(i in range(length(input_gds_files))){
        String output_vcf = basename(output_vcf_files[i], ".gz")
        # Convert files
        call SEQARRAY.convert_seq_gds_to_vcf as convert_seq_gds_to_vcf {
            input:
                in_seq_gds = input_gds_files[i],
                out_vcf = output_vcf,
                seq_gds_variant_annot = input_gds_variant_annot[i],
                variant_ids = variant_ids,
                variant_annot_xref_col = variant_annot_xref_col,
                variant_annot_gds_id_col = variant_annot_gds_id_col,
                seq_gds_sample_annot = input_gds_sample_annot[i],
                sample_ids = sample_ids,
                sample_annot_xref_col = sample_annot_xref_col,
                sample_annot_gds_id_col = sample_annot_gds_id_col,
                cpu = cpu,
                mem_gb = mem_gb
        }
    }

    output{
        Array[File] results_files = convert_seq_gds_to_vcf.output_vcf
    }
}

