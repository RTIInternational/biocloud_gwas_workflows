version 1.1

import "pulmonary_longitudinal_gxe.wdl" as GXE
import "rti-tsv-utils.wdl" as TSV

workflow pulmonary_longitudinal_gxe_wf{

    input {

        Array[File] file_bgens
        Array[File] file_bgis
        File file_sample
        File file_pheno
        String ancestry
        String pheno
        String omega3
        Int? chunk_size
        String file_out_prefix
        
        # Resources
        Int cpu = 2
        Int mem_gb = 8
        String image_source = "docker"
        String? ecr_repo

    }

    scatter(i in range(length(file_bgens))) {

        call GXE.run_gxe {
            input:
                file_bgen = file_bgens[i],
                file_bgi = file_bgis[i],
                file_sample = file_sample,
                file_pheno = file_pheno,
                ancestry = ancestry,
                pheno = pheno,
                omega3 = omega3,
                file_out_prefix = "~{file_out_prefix}_~{i}",
                cpu = cpu,
                mem_gb = mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

    }

    call TSV.tsv_append as cat_gxe_results{
        input:
            input_files = run_gxe.gxe_results,
            output_prefix = "~{file_out_prefix}",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        Array[File] gxe_results = cat_gxe_results.out_tsv
    }

}
