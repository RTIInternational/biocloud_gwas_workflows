import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT
import "biocloud_gwas_workflows/biocloud_wdl_tools/rti-tsv-utils/rti-tsv-utils.wdl" as SORT

workflow convert_variant_ids_bfile_wf{
    File bed_in
    File bim_in
    File fam_in
    Array[File] ref_files
    Array[String] chrs
    String output_basename

    String? in_sep
    String? in_missing_allele
    String? in_deletion_allele
    String? ref_deletion_allele
    Int? in_chunk_size
    Int? ref_chunk_size
    Boolean? rescue_rsids
    String? output_compression

    # PLINK resources
    Int plink_cpu = 1
    Int plink_mem_gb = 2

    # Convert resources
    Int convert_cpu = 1
    Int convert_mem_gb = 4

    String container_source = "docker"

    # Parallelize
    scatter(i in range(length(chrs))){

        # Split bfile
        call PLINK.make_bed as split_bfile{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in= fam_in,
                chr = "${chrs[i]}",
                output_basename = output_basename,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb,
                container_source = container_source
        }

        # Convert IDs
        call IDCONVERT.convert_variant_ids{
            input:
                chr = "${chrs[i]}",
                in_file = split_bfile.bim_out,
                in_header = 0,
                in_sep = in_sep,
                in_id_col = 1,
                in_chr_col = 0,
                in_pos_col = 3,
                in_a1_col = 4,
                in_a2_col = 5,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                in_chunk_size = in_chunk_size,
                ref = ref_files[i],
                ref_deletion_allele = ref_deletion_allele,
                ref_chunk_size = ref_chunk_size,
                output_filename = "${output_basename}_chr${chrs[i]}.bim",
                rescue_rsids = rescue_rsids,
                output_compression = output_compression,
                cpu = convert_cpu,
                mem_gb = convert_mem_gb,
                container_source = container_source
        }

    }

    call PLINK.merge_beds as merge_chrs{
        input:
            bed_in = split_bfile.bed_out,
            bim_in = convert_variant_ids.output_file,
            fam_in = split_bfile.fam_out,
            output_basename = output_basename,
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            container_source = container_source
    }

    output{
        File bed_out = merge_chrs.bed_out
        File bim_out = merge_chrs.bim_out
        File fam_out = merge_chrs.fam_out
        File log_out = merge_chrs.plink_log
    }
}
