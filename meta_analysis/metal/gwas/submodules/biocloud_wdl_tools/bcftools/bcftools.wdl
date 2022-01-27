task view{
    # Utility for splitting a VCF file into chunks of N variants
    File vcf_in
    File? vcf_tbi_in

    File? samples_file
    String? regions
    File? regions_file
    String? targets
    File? targets_file
    Float? maf_filter
    Int? min_ac

    String output_filename
    String output_type

    # Runtime environment
    String docker = "rtibiocloud/bcftools:v1.9-8875c1e"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        # Merge into single VCF
        bcftools view \
            --threads ${cpu} \
            -o ${output_filename} \
            -O ${output_type} \
            ${'-S ' + samples_file + " --force-sample"} \
            ${'-r ' + regions} \
            ${'-R ' + regions_file} \
            ${'-t ' + targets} \
            ${'-T ' + targets_file} \
            ${'-q ' + maf_filter + ":minor"} \
            ${'--min-ac ' + min_ac + ":minor"} \
            ${vcf_in}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File vcf_out = "${output_filename}"
    }
}
