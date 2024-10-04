version 1.1

struct FULL_QC_GENOTYPE_FILES {
    Array[File] beds
    Array[File] bims
    Array[File] fams
}

struct BASIC_QC_GENOTYPE_FILES {
    Array[File] beds
    Array[File] bims
    Array[File] fams
}

struct VARIANT_LISTS {
    File initial
    File final
    File low_call_rate
    File failed_hwe
}

struct SAMPLE_LISTS {
    File initial
    File final_full_qc
    File final_basic_qc
    File low_call_rate
    File excess_homozygosity
    File related
    File failed_sex_check
}

struct COUNTS {
    Int variant_initial
    Int variant_final
    Int variant_low_call_rate
    Int variant_failed_hwe
    Int sample_initial
    Int sample_final_full_qc
    Int sample_final_basic_qc
    Int sample_low_call_rate
    Int sample_excess_homozygosity
    Int sample_related
    Int sample_failed_sex_check
}

struct RELATEDNESS_WF_OUTPUTS {
    File annotated_kinships
    File kinship_id_maps
}

struct SEX_CHECK_WF_OUTPUTS {
    File sex_check_report
}
