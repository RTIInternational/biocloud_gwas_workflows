import "biocloud_gwas_workflows/biocloud_wdl_tools/split_utils/split_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/rvtests/rvtests.wdl" as RV
import "biocloud_gwas_workflows/biocloud_wdl_tools/make_gwas_summary_stats/make_gwas_summary_stats.wdl" as STAT
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV

task strip_rvtests_headers{
    # Utility module to remove comment lines from header of RVTests files
    # Special module needed so we can join RVtests from parallel chunks (tsv-utils doesn't accept tsvs with comments before header)
    File sumstats_in
    String output_basename
    String output_file = output_basename + ".header_stripped.tsv"

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command<<<
        input_file=${sumstats_in}

        # Unzip tsv input file if necessary
        if [[ ${sumstats_in} =~ \.gz$ ]]; then
            log_info "${sumstats_in} is gzipped. Unzipping..."
            gunzip -c ${sumstats_in} > input.txt
            input_file=input.txt
        fi

        grep -v "^#" $input_file > ${output_file}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File  sumstats_out = "${output_file}"
    }
}

workflow rvtests_chr_wf{
    File vcf_in
    File info_in
    File pheno_file
    String pheno_name
    File covar_file
    Array[String] covars
    File? kinship_matrix
    String dosage

    # For annotating with population MAF info
    File pop_maf_file
    String maf_population

    # Optional for X chr
    File? xHemiKinship_matrix
    Boolean? xHemi

    # Number of records per VCF split
    Int records_per_split = 50000

    # Set output basename
    String? user_output_basename
    String default_output_basename = basename(sub(vcf_in, "\\.gz$",""), ".vcf")
    String output_basename = select_first([user_output_basename, default_output_basename])

    # Chunk VCF and INFO files for parallel processing
    call SPLIT.split_vcf as split_vcf{
        input:
            input_vcf = vcf_in,
            records_per_split = records_per_split,
            output_basename = output_basename
    }

    # Loop through splits and do association testing on each
    scatter(split_index in range(length(split_vcf.split_vcfs))){

        String split_output_basename = basename(sub(split_vcf.split_vcfs[split_index], "\\.gz$",""), ".vcf")

        # Run rvtests for association
        call RV.rvtests{
            input:
                inVCF = split_vcf.split_vcfs[split_index],
                phenoFile = pheno_file,
                output_basename = split_output_basename,
                covarFile = covar_file,
                phenoName = pheno_name,
                covarsMaybe = covars,
                kinship = kinship_matrix,
                xHemiKinship = xHemiKinship_matrix,
                xHemi = xHemi
        }

        # Remove header from association output file
        File assoc_file = rvtests.assoc_outputs[0]
        call strip_rvtests_headers{
            input:
                sumstats_in = assoc_file,
                output_basename = basename(assoc_file, ".gz")
        }
    }

    # Collect chunked sumstats files into single zip folder
    call COLLECT.collect_large_file_list_wf as collect_sumstats{
        input:
      #      input_files = flatten_sumstats.flat_array,
            input_files = strip_rvtests_headers.sumstats_out,
            output_dir_name = output_basename + ".rvtests_output"
    }

    # Concat all sumstats files into single sumstat file
    call TSV.tsv_append as cat_sumstats{
        input:
            tsv_inputs_tarball = collect_sumstats.output_dir,
            output_filename = output_basename + ".rvtests.merged_assoc.tsv"
    }

    call STAT.make_gwas_summary_stats as annotate_sumstats{
        input:
            file_in_summary_stats = cat_sumstats.tsv_output,
            file_in_info = info_in,
            file_in_pop_mafs = pop_maf_file,
            file_in_summary_stats_format = "rvtests",
            population = maf_population,
            file_out_prefix = output_basename + ".formatted"
    }

    output{
        File summary_stats = annotate_sumstats.output_file
    }
}