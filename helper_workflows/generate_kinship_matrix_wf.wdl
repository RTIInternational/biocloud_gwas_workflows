import "biocloud_gwas_workflows/biocloud_wdl_tools/merge_utils/merge_utils.wdl" as MERGE
import "biocloud_gwas_workflows/biocloud_wdl_tools/rvtests/rvtests.wdl" as RV
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

workflow generate_kinship_matrix_wf{
    Array[File] input_vcfs
    Array[String] chrs
    File ped_file
    String output_basename

    # RVTests options
    Boolean xHemi
    Boolean useBaldingNicols
    Boolean useIBS
    String? dosage
    String? xLabel
    Float? maxMiss
    Float? minMAF
    Float? minSiteQual

    Int split_kinship_cpu = 6
    Int split_kinship_mem_gb = 6

    Int combine_kinship_cpu = 16
    Int combine_kinship_mem_gb = 16

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

    # Remove empty kinship files
    # WDL hack because WDL doesn't support optional output and both vcf2kinship output files are technically optional
    # e.g. Kinship files won't be produced for X chr; xHemiKinship files won't exist for autosomes
    # To get around this, vcf2kinship 'touches' both output files so they technically exist but in some cases are just empty placeholders
    # These files need to be removed
    call UTILS.remove_empty_files as remove_empty_auto{
        input:
            input_files = make_split_kinship.kinship_matrix
    }

    call UTILS.remove_empty_files as remove_empty_x{
        input:
            input_files = make_split_kinship.xHemi_kinship_matrix
    }

    # Merge scattered kinship autosomal matrices
    call RV.combineKinship as combine_kin{
        input:
            kinship_matrices = remove_empty_auto.non_empty_files,
            output_basename = "${output_basename}.merged.split.auto",
            cpu = combine_kinship_cpu,
            mem_gb = combine_kinship_mem_gb
    }

    # This is kind of a stupid trick to try and get
    File? null_xKin
    File? xHemiKinMat = if(length(remove_empty_x.non_empty_files) > 0) then remove_empty_x.non_empty_files[0] else null_xKin

    output{
        #File kinship_matrix = combine_kin.kinship_matrix
        #File? xHemi_kinship_matrix = select_first(make_split_kinship.xHemi_kinship_matrix)
        File kinship_matrix = combine_kin.kinship_matrix
        File? xHemi_kinship_matrix = xHemiKinMat
    }
}