import "biocloud_gwas_workflows/biocloud_wdl_tools/merge_utils/merge_utils.wdl" as MERGE
import "biocloud_gwas_workflows/biocloud_wdl_tools/rvtests/rvtests.wdl" as RV

workflow generate_kinship_matrix_wf{
    Array[File] input_vcfs
    Array[String] chrs
    File ped_file
    String output_basename

    # Merge VCF options
    Int merge_vcf_cpus = 8
    Int merge_vcf_mem_gb = 16

    # RVTests options
    Boolean xHemi
    Boolean useBaldingNicols
    Boolean useIBS
    String? dosage
    String? xLabel
    Float? maxMiss
    Float? minMAF
    Float? minSiteQual

    Int giant_kinship_cpu = 32
    Int giant_kinship_mem_gb = 32

    Int split_kinship_cpu = 6
    Int split_kinship_mem_gb = 12

    Int combine_kinship_cpu = 16
    Int combine_kinship_mem_gb = 16

    # Merge VCFs into single large VCF
    call MERGE.merge_vcfs{
        input:
            input_vcfs = input_vcfs,
            output_filename = "${output_basename}.merged.vcf.gz"
    }

    # Generate kinship matrix from merged VCF
    call RV.vcf2kinship as make_giant_kinship{
        input:
            inputVcf = merge_vcfs.merged_vcf,
            dosage = dosage,
            xHemi = xHemi,
            xLabel = xLabel,
            maxMiss = maxMiss,
            minMAF = minMAF,
            minSiteQual = minSiteQual,
            useBaldingNicols = useBaldingNicols,
            useIBS = useIBS,
            output_basename = "${output_basename}.giant",
            cpu = giant_kinship_cpu,
            mem_gb = giant_kinship_mem_gb
    }

    # Do scattered workflow
    scatter(chr_index in range(length(input_vcfs))){

        # Only set xHemi and pass pedfile for X chromosome
        String chr = chrs[chr_index]
        Boolean xHemi_chr = if((chr == "X") && (xHemi)) then true else false
        File? null_pedfile
        File? pedfile_chr = if((chr == "X") && (xHemi)) then ped_file else null_pedfile

        call RV.vcf2kinship as make_split_kinship{
            input:
                inputVcf = input_vcfs[chr_index],
                pedfile = pedfile_chr,
                dosage = dosage,
                xHemi = xHemi_chr,
                xLabel = xLabel,
                maxMiss = maxMiss,
                minMAF = minMAF,
                minSiteQual = minSiteQual,
                useBaldingNicols = useBaldingNicols,
                useIBS = useIBS,
                output_basename = "${output_basename}.split.${chr}",
                cpu = split_kinship_cpu,
                mem_gb = split_kinship_mem_gb
        }
    }

    # Merge scattered kinship autosomal matrices
    call RV.combineKinship as combine_kin{
        input:
            kinship_matrices = make_split_kinship.kinship_matrix,
            output_basename = "${output_basename}.merged.split.auto",
            cpu = combine_kinship_cpu,
            mem_gb = combine_kinship_mem_gb
    }

    output{
        File giant_kinship_matrix = make_giant_kinship.kinship_matrix
        File? giant_xHemi_kinship_matrix = make_giant_kinship.xHemi_kinship_matrix
        File split_kinship_matrix = combine_kin.kinship_matrix
        File? split_xHemi_kinship_matrix = select_first(make_split_kinship.xHemi_kinship_matrix)
    }
}