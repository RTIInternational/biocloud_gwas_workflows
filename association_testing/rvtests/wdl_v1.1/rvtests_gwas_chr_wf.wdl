version 1.1

import "split_utils.wdl" as SPLIT
import "rvtests.wdl" as RV
import "make_gwas_summary_stats.wdl" as STAT
import "rti-tsv-utils.wdl" as TSV

workflow rvtests_gwas_chr_wf{

    input {

        # Input VCF and INFO files
        File vcf_in
        File info_in
        String info_format = "info"
        File pheno_file
        String pheno_name
        File? covar_file
        Array[String]? covars

        String? dosage

        # For annotating with population MAF info
        File? pop_maf_file

        # Optional for X chr
        File? kinship
        File? xHemiKinship
        Boolean? xHemi

        Array[String]? singleTestsMaybe
        Array[String]? burdenTestsMaybe
        Array[String]? vtTestsMaybe
        Array[String]? kernelTestsMaybe
        Array[String]? metaTestsMaybe

        Boolean? inverseNormal
        Boolean? useResidualAsPhenotype
        Boolean? sex
        Boolean? qtl
        Boolean? multipleAllele
        String? xLabel
        Int rvtests_cpu = 1
        Int rvtests_mem_gb = 2

        # Number of records per VCF split
        Int records_per_split = 50000
        Int split_vcf_cpu = 8
        Boolean split_vcf_records = true

        # Set output basename
        String? user_output_basename
        String default_output_basename = basename(sub(vcf_in, "\\.gz$",""), ".vcf")
        String output_basename = select_first([user_output_basename, default_output_basename])

        # Container settings
        String image_source = "docker"
        String? ecr_repo

    }

    # Optionally chunk VCF and INFO files for parallel processing
    if(split_vcf_records){
        call SPLIT.split_vcf as split_vcf{
            input:
                input_vcf = vcf_in,
                records_per_split = records_per_split,
                output_basename = output_basename,
                cpu = split_vcf_cpu,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Loop through splits and do association testing on each
    Array[File] split_vcfs = select_first([split_vcf.split_vcfs, [vcf_in]])

    scatter(split_index in range(length(split_vcfs))){

        String split_output_basename = basename(sub(split_vcfs[split_index], "\\.gz$",""), ".vcf")

        # Run rvtests for association
        call RV.rvtests{
            input:
                inVCF = split_vcfs[split_index],
                phenoFile = pheno_file,
                output_basename = split_output_basename,
                covarFile = covar_file,
                phenoName = pheno_name,
                covarsMaybe = covars,
                kinship = kinship,
                xHemiKinship = xHemiKinship,
                xHemi = xHemi,
                singleTestsMaybe = singleTestsMaybe,
                burdenTestsMaybe = burdenTestsMaybe,
                vtTestsMaybe = vtTestsMaybe,
                kernelTestsMaybe = kernelTestsMaybe,
                metaTestsMaybe = metaTestsMaybe,
                inverseNormal = inverseNormal,
                useResidualAsPhenotype = useResidualAsPhenotype,
                sex = sex,
                qtl = qtl,
                multipleAllele = multipleAllele,
                xLabel = xLabel,
                dosage = dosage,
                cpu = rvtests_cpu,
                mem_gb = rvtests_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Remove header from association output file
        File assoc_file = rvtests.assoc_outputs[0]
        call strip_rvtests_headers{
            input:
                sumstats_in = assoc_file,
                output_basename = basename(assoc_file, ".gz"),
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Gather chunked RVTests output if > 1 split
    if(length(split_vcfs) > 1){
        # Concat all sumstats files into single sumstat file
        call TSV.tsv_append as cat_rvtests_sumstats{
            input:
                input_files = strip_rvtests_headers.sumstats_out,
                output_prefix = output_basename + ".rvtests.MetaAssoc",
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    # Use the combined sumstats file if >1 splits was merged, otherwise just use output from strip_headers call
    File full_sumstats = select_first([cat_rvtests_sumstats.out_tsv, strip_rvtests_headers.sumstats_out[0]])

    # Annotate sumstats with features from info file (R2, MAF) and pop MAF file (MAF from pop of interest)
    call STAT.make_gwas_summary_stats as annotate_sumstats{
        input:
            file_in_summary_stats = full_sumstats,
            file_in_info = info_in,
            file_in_pop_mafs = pop_maf_file,
            file_in_summary_stats_format = "rvtests",
            file_in_info_format = info_format,
            file_out_prefix = output_basename + ".formatted",
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    output{
        File summary_stats = annotate_sumstats.output_file
    }
}

task strip_rvtests_headers{
    # Utility module to remove comment lines from header of RVTests files
    # Special module needed so we can join RVtests from parallel chunks (tsv-utils doesn't accept tsvs with comments before header)

    input {

        File sumstats_in
        String output_basename
        String output_file = output_basename + ".header_stripped.tsv"

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:a6d2b38300ce017add71440577d5b0a90460d0e57fd7aec21dd0d1b0761bbfb2"
        String ecr_image = "ubuntu:22.04_19478ce7fc2ff"
        String image_source = "docker"
        String? ecr_repo
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 1

    }

    command<<<
        set -e
        input_file=~{sumstats_in}

        # Unzip tsv input file if necessary
        if [[ ~{sumstats_in} =~ \.gz$ ]]; then
            echo "~{sumstats_in} is gzipped. Unzipping..."
            gunzip -c ~{sumstats_in} > input.txt
            input_file=input.txt
        fi

        grep -v "^#" $input_file > ~{output_file}
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File  sumstats_out = "~{output_file}"
    }
}
