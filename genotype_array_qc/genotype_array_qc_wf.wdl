import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/impute2_id_conversion/impute2_id_conversion_wf.wdl" as IDCONVERT

workflow genotype_array_qc_wf{
    File bed
    File bim
    File fam
    Array[String] chrs
    Array[File] id_legend_files
    String output_basename

    # PAR/NONPAR Split/Merge parameters
    String build_code
    Boolean no_fail = true

    # Sample call-rate filter cutoffs
    Float failed_sample_missed_call_rate = 0.99
    Float min_sample_call_rate

    # Sex-check filter parameters
    Boolean filter_discrepant_sex
    Float max_female_f = 0.2
    Float min_male_f = 0.8

    # Kinship filter parameters
    Boolean filter_related_samples
    Float max_relatedness_coefficient

    # Admixed ancestry filtering parameters
    Boolean filter_admixed_samples
    Float max_pop_admixture_threshold

    # Variant filtering cutoffs
    Float min_site_call_rate
    Float hwe_filter_pvalue
    Float max_homozygosity

    # Resources for various modules
    Int split_bed_cpu = 4
    Int split_bed_mem_gb = 8
    Int merge_bed_cpu = 2
    Int merge_bed_mem_gb = 12

    # Remove phenotype from fam file
    call PLINK.remove_fam_phenotype{
        input:
            fam_in = fam,
            output_basename = "${output_basename}.no_pheno"
    }

    # Convert variant IDs to impute2 format and remove duplicate variants
    call IDCONVERT.impute2_id_conversion_wf as convert_impute2_ids{
        input:
            bed_in = bed,
            bim_in = bim,
            fam_in = remove_fam_phenotype.fam_out,
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
            duplicate_id_mem_gb = split_bed_mem_gb
    }

    # Remove failed subjects with >99% missing data
    call PLINK.make_bed as filter_failed_samples{
        input:
            bed_in = convert_impute2_ids.bed_out,
            bim_in = convert_impute2_ids.bim_out,
            fam_in = convert_impute2_ids.fam_out,
            output_basename = "${output_basename}.filter_failed_samples",
            mind = failed_sample_missed_call_rate
    }

    # TODO: STRUCTURE WF to determine ancestry
    # TODO: Partition data by ancestry

    # For each ancestry group, do the following filtering steps
    # TODO: SNP Call rate filter
    # TODO: HWE filter
    # TODO: Set het haploids to missing
    # TODO: Subject call rate filter (autosomes)
    # TODO: Excess homozygosity filtering
    # TODO: Relatedness workflow
    # TODO: Sex check WF
    # TODO: Optionally Remove samples based on relatedness
    # TODO: Optionally Remove samples based on discrepant sex

    # TODO: Merge PAR/NONPAR regoins of chrX
    # TODO: Flag individuals missing chrX or other chr

    output{
        File bed_out = merge_beds.bed_out
        File bim_out = merge_beds.bim_out
        File fam_out = merge_beds.fam_out
        File fam_out_no_pheno = remove_fam_phenotype.fam_out
    }
}