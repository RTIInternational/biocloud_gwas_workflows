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

    File? bed
    File? bim
    File? fam

    String output_basename
    Float maf = 0.01

    String build_code
    Array[String] chrs
    Array[File] id_legend_files
    Boolean remove_duplicates = true
    String in_monomorphic_allele = "0"
    String in_deletion_allele = "-"
    String ref_deletion_allele = "."

    # Resources for various modules
    Int convert_bed_cpu = 6
    Int convert_bed_mem_gb = 16

    Int plink_cpu = 8
    Int plink_mem_gb = 16


    # Convert pgen/pvar file to bed/bim/fam format
    if(!defined(bed)){
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
    }

    File converted_bed = select_first([bed, convert_to_bed.bed_out])
    File converted_bim = select_first([bim, convert_to_bed.bim_out])
    File converted_fam = select_first([fam, convert_to_bed.fam_out])

    # Convert variant IDs to impute2 format and remove duplicate variants
    call IDCONVERT.impute2_id_conversion_wf as convert_impute2_ids{
        input:
            bed_in = converted_bed,
            bim_in = converted_bim,
            fam_in = converted_fam,
            output_basename = output_basename,
            chrs = chrs,
            id_legend_files = id_legend_files,
            remove_duplicates = true,
            build_code = build_code,
            plink_cpu = plink_cpu,
            plink_mem_gb = plink_mem_gb,
            id_convert_cpu = 1,
            id_convert_mem_gb = 6,
            in_monomorphic_allele = in_monomorphic_allele,
            in_deletion_allele = in_deletion_allele,
            ref_deletion_allele = "."
    }

    output{
        File bed_out = convert_impute2_ids.bed_out
        File bim_out = convert_impute2_ids.bim_out
        File fam_out = convert_impute2_ids.fam_out
    }
}

