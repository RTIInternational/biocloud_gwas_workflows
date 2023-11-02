version 1.1

import "tasks/liftover_utils.wdl" as utils

workflow liftOver {

meta {
    description: "Transfer genomic coordinates from one genome build to another (liftOver) using the liftover function from the MungeSumstats package by Neurogenomics Lab."
    maintainer: "Jesse Marks <jmarks@rti.org>"
    software_version: "1.7.10"
    software_website: "https://neurogenomics.github.io/MungeSumstats/reference/liftover.html"
    }

    input {
        File file_name
        String output_name
        String sep
        String snp_name
        String chrom_name
        String pos_name
        String ref_genome
        String convert_ref_genome
        String chain_source

        String docker = "rtibiocloud/neurogenomicslab_mungesumstats:1.7.10_3d50aed"
        Int cpu = 1
        Int mem = ceil(size(file_name, "GB"))
    }

    call utils.liftOverTask as lift {
        input:
            file_name = file_name,
            output_name = output_name,
            sep = sep,
            snp_name = snp_name,
            chrom_name = chrom_name,
            pos_name = pos_name,
            ref_genome = ref_genome,
            convert_ref_genome = convert_ref_genome,
            chain_source = chain_source,
            docker = docker,
            cpu = cpu,
            mem = mem
    }

    output {
        File final_file = lift.output_file
    }
}
