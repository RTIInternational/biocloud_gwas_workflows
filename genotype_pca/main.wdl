import "biocloud_gwas_workflows/genotype_pca/utils.wdl" as UTILS

workflow eigenstrat_smartpca{

    String study_name
    String ancestry

    File bedfile
    File bimfile
    File famfile

    String docker_ubuntu = "ubuntu:22.04"
    String docker_plink1_9 = "rtibiocloud/plink:v1.9_178bb91"
    String docker_eigensoft = "rtibiocloud/eigensoft:v6.1.4_2d0f99b"
    String cpu = 1
    String mem = 2

    call UTILS.locate_high_ld_regions {
        input:
            bimfile = bimfile,
            docker = docker_ubuntu,
            cpu = cpu,
            mem = mem
    }

    call UTILS.remove_high_ld_regions {
        input:
            study_name = study_name,
            ancestry = ancestry,
            bedfile = bedfile, 
            bimfile = bimfile,
            famfile = famfile, 
            high_ld_regions = locate_high_ld_regions.high_ld_regions,
            docker = docker_plink1_9,
            cpu = cpu,
            mem = mem
    }


    call UTILS.ld_pruning {
        input:
            study_name = study_name,
            ancestry = ancestry,
            bedfile = remove_high_ld_regions.output_bed, 
            bimfile = remove_high_ld_regions.output_bim,
            famfile = remove_high_ld_regions.output_fam, 
            docker = docker_plink1_9,
            cpu = cpu,
            mem = mem
    }

    call UTILS.merge_pruned {
        input:
            pruned_files = ld_pruning.pruned_files,
            docker = docker_ubuntu,
            cpu = cpu,
            mem = mem
    }

    call UTILS.extract_ld_variants {
        input:
            study_name = study_name,
            ancestry = ancestry,
            bedfile = remove_high_ld_regions.output_bed,
            bimfile = remove_high_ld_regions.output_bim,
            famfile = remove_high_ld_regions.output_fam,
            combined_variants = merge_pruned.combined_variants,
            docker = docker_plink1_9,
            cpu = cpu,
            mem = mem
    }

    call UTILS.rename_bimfam {
        input:
            bimfile = extract_ld_variants.output_bim,
            famfile = extract_ld_variants.output_fam,
            docker = docker_plink1_9,
            cpu = cpu,
            mem = mem
    }

    call UTILS.run_smartpca {
        input:
            ancestry = ancestry,
            study_name = study_name,
            bedfile = extract_ld_variants.output_bed,
            bimfile = rename_bimfam.ids_renamed_bim,
            famfile = rename_bimfam.ids_renamed_fam,
            docker = docker_eigensoft,
            cpu = cpu,
            mem = mem
    }

    call UTILS.create_final_file {
        input:
            study_name = study_name,
            ancestry = ancestry,
            evec_file = run_smartpca.eigenvectors,
            famfile = famfile,

            docker = docker_ubuntu,
            cpu = cpu,
            mem = mem
    }

    output {
        File final_file = create_final_file.final_file
    }
}
