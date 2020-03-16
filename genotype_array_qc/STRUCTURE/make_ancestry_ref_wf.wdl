import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/impute2_id_conversion/impute2_id_conversion_wf.wdl" as IDCONVERT

task convert_to_bed{
    File pgen_in
    File pvar_in
    File psam_in
    String output_basename
    String input_prefix = basename(pgen_in, ".pgen.zst")
    Float maf = 0.01

    String docker = "rtibiocloud/plink:v2.0-8875c1e"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command{
        set -e
        plink2 --zst-decompress ${pgen_in} > ${input_prefix}.pgen
        plink2 --zst-decompress ${pvar_in} > ${input_prefix}.pvar
        ln -s ${psam_in} ${input_prefix}.psam

        plink2 --pfile ${input_prefix} \
            --make-bed \
            --threads ${cpu} \
            --snps-only just-acgt \
            --max-alleles 2 \
            --autosome \
            --maf ${maf} \
            --out ${output_basename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File bed_out = "${output_basename}.bed"
        File bim_out = "${output_basename}.bim"
        File fam_out = "${output_basename}.fam"
        File plink_log = "${output_basename}.log"
    }
}

workflow make_ancestry_ref_wf{
    File pgen_zst
    File pvar_zst
    File psam_zst

    Array[String] chrs
    Array[File] id_legend_files

    String output_basename
    String build_code
    Float maf = 0.01

    # Resources for various modules
    Int convert_bed_cpu = 6
    Int convert_bed_mem_gb = 16
    Int split_bed_cpu = 4
    Int split_bed_mem_gb = 8
    Int merge_bed_cpu = 8
    Int merge_bed_mem_gb = 16
    Int x_merge_cpu = 8
    Int x_merge_mem_gb = 16


    # Convert pgen/pvar file to bed/bim/fam format
    call convert_to_bed{
        input:
            pgen_in = pgen_zst,
            pvar_in = pvar_zst,
            psam_in = psam_zst,
            maf = maf,
            cpu = convert_bed_cpu,
            mem_gb = convert_bed_mem_gb,
            output_basename = "${output_basename}.raw"
    }

    # Convert variant IDs to impute2 format and remove duplicate variants
    call IDCONVERT.impute2_id_conversion_wf as convert_impute2_ids{
        input:
            bed_in = convert_to_bed.bed_out,
            bim_in = convert_to_bed.bim_out,
            fam_in = convert_to_bed.fam_out,
            output_basename = output_basename,
            chrs = chrs,
            id_legend_files = id_legend_files,
            remove_duplicates = true,
            build_code = build_code,
            no_fail = true,
            split_bed_cpu = split_bed_cpu,
            split_bed_mem_gb = split_bed_mem_gb,
            merge_bed_cpu = merge_bed_cpu,
            merge_bed_mem_gb = merge_bed_mem_gb,
            duplicate_id_cpu = split_bed_cpu,
            duplicate_id_mem_gb = split_bed_mem_gb,
            x_merge_cpu = x_merge_cpu,
            x_merge_mem_gb = x_merge_mem_gb
    }
    output{
        File bed_out = convert_impute2_ids.bed_out
        File bim_out = convert_impute2_ids.bim_out
        File fam_out = convert_impute2_ids.fam_out
    }
}

