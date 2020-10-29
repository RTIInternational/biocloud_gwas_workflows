import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

task merge_samples_to_remove{
    # Return whether input chrs are in numerical sort order
    Array[File] input_files
    String output_filename
    command<<<
        cat ${sep=" " input_files} | sort | uniq > ${output_filename}
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "100 MB"
    }

    output{
        File output_file = "${output_filename}"
    }
}

workflow missing_chr_filter_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    # HWE p-value filter threshold
    Float max_sample_chr_missing_rate = 0.99

    # Workflow can optionally be used on a single chr
    Array[String]? chrs

    Int plink_chr_cpu = 1
    Int plink_chr_mem_gb = 2
    Int plink_filter_cpu = 1
    Int plink_filter_mem_gb = 2

    # Get all chromosomes in bim if no list provided by user
    if(!defined(chrs)){
        call PLINK.get_bim_chrs{
            input:
                bim_in = bim_in
        }
    }

    Array[String] scatter_chrs = select_first([chrs, get_bim_chrs.chrs])

    scatter(chr in scatter_chrs){

        # Identify any samples that are missing most genotypes for a chr
        call PLINK.get_samples_missing_chr as filter_missing{
            input:
                bed_in = bed_in,
                bim_in = bim_in,
                fam_in = fam_in,
                chr = chr,
                missing_threshold = max_sample_chr_missing_rate,
                output_basename = output_basename,
                cpu = plink_chr_cpu,
                mem_gb = plink_chr_mem_gb
        }
    }

    # Combine snps to remove across all chrs
    call merge_samples_to_remove{
        input:
            input_files = filter_missing.missing_samples,
            output_filename = "${output_basename}.missing.merged.remove"
    }

    # Remove failed samples for any chr from original dataset
    call PLINK.make_bed as filter_samples_missing_chr{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            remove_samples = merge_samples_to_remove.output_file,
            chrs = scatter_chrs,
            output_basename = output_basename,
            cpu = plink_filter_cpu,
            mem_gb = plink_filter_mem_gb
    }

    output{
        File bed_out = filter_samples_missing_chr.bed_out
        File bim_out = filter_samples_missing_chr.bim_out
        File fam_out = filter_samples_missing_chr.fam_out
        File samples_missing_chr = merge_samples_to_remove.output_file
    }
}
