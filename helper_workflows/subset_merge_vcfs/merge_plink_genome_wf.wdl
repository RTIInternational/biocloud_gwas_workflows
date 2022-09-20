import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow merge_plink_wf{
    Array[File] beds_in
    Array[File] bims_in
    Array[File] fams_in

    String output_basename
    Int cpu
    Int mem_gb

    call PLINK.merge_beds{
        input:
            bed_in = beds_in,
            bim_in = bims_in,
            fam_in = fams_in,
            output_basename = output_basename,
            cpu = cpu,
            mem_gb = mem_gb
    }
    
    output {
        File bed = merge_plink.bed_out
        File bim = merge_plink.bim_out
        File fam = merge_plink.fam_out
        File log = merge_plink.plink_log
    }
    
}

