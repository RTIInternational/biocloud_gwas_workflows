version 1.1

import "merge_utils.wdl" as MERGE
import "rvtests.wdl" as RV
import "utils.wdl" as UTILS

workflow generate_kinship_matrix_wf{

    input {
        Array[File] input_vcfs
        Array[String] chrs
        File ped_file
        String output_basename

        # RVTests options
        Boolean useBaldingNicols
        Boolean useIBS
        String? dosage
        Float? maxMiss
        Float? minMAF
        Float? minSiteQual

        # Sex chromosome options
        Boolean xHemi = false
        String xLabel = "X"

        # Runtime/Resource options
        Int split_kinship_cpu = 1
        Int split_kinship_mem_gb = 2
        Int combine_kinship_cpu = 2
        Int combine_kinship_mem_gb = 4

        # Placeholder null pedfile because we only want to pass the pedfile when we're doing xHemi
        # Really this is a result of a poorly designed WDL module
        File? null_pedfile

        # Container
        String image_source = "docker"
        String? ecr_repo
    }

    # Generate kinship matrix for each chromosome in parallel
    scatter(chr_index in range(length(input_vcfs))){

        # Only set xHemi and pass pedfile for X chromosome
        String chr = chrs[chr_index]
        Boolean xHemi_chr = if((chr == xLabel) && (xHemi)) then true else false
        File? pedfile_chr = if((chr == xLabel) && (xHemi)) then ped_file else null_pedfile

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
                output_basename = "~{output_basename}.split.chr.~{chr}",
                cpu = split_kinship_cpu,
                mem_gb = split_kinship_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Remove empty kinship files
    # WDL hack because WDL doesn't support optional output and both vcf2kinship output files are technically optional
    # e.g. Kinship files won't be produced for X chr; xHemiKinship files won't exist for autosomes
    # To get around this, vcf2kinship 'touches' both output files so they technically exist but in some cases are just empty placeholders
    # These files need to be removed
    call UTILS.remove_empty_files as remove_empty_auto{
        input:
            input_files = make_split_kinship.kinship_matrix,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Only try to get the xHemi kinship matrix if xHemi workflow input is true
    if(xHemi){
        call UTILS.remove_empty_files as remove_empty_x{
            input:
                input_files = make_split_kinship.xHemi_kinship_matrix,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
        File xHemiKinshipMat = remove_empty_x.non_empty_files[0]
    }

    # Merge scattered kinship autosomal matrices
    call RV.combineKinship as combine_kin{
        input:
            kinship_matrices = remove_empty_auto.non_empty_files,
            vcf2kinship_logs = make_split_kinship.kinship_log,
            output_basename = "~{output_basename}.merged.split.auto",
            cpu = combine_kinship_cpu,
            mem_gb = combine_kinship_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }


    output{
        File kinship = combine_kin.kinship_matrix
        File? xHemiKinship = xHemiKinshipMat
    }
}