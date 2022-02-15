import "biocloud_gwas_workflows/liftover_genomic_annotations/liftover_utils.wdl" as UTILS


workflow genome_liftover {
    String chainfile =  "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
    File input_sumstats
    Int chromosome_col  # zero-based column index
    Int position_col    # zero-based column index
    String final_file

    String docker = "rtibiocloud/liftover:v423_1471834"  

    call UTILS.download_chainfile as chain {
        input:
            chainfile = chainfile,
            docker = docker
    }

    call UTILS.create_bedfile as bed {
        input:
            input_sumstats = input_sumstats,
            chromosome_col = chromosome_col,
            position_col = position_col,
            docker = docker

    }

    call UTILS.perform_liftover as lift {
        input:
            input_bedfile = bed.sumstats_bed,
            chainfile = chain.unzipped_output,
            docker = docker
    }

    call UTILS.final_sumstats  as final {
        input:
            original_sumstats = input_sumstats,
            new_bed = lift.output_bed,
            unmapped_bed = lift.unmapped_bed,
            output_name = final_file,
            position_col = position_col,
            chr_col = chromosome_col,
            docker = docker
    }

    output {
        File mapped_stats = final.output_sumstats
        File unmapped_stats = final.unmapped_sumstats
    }
}
