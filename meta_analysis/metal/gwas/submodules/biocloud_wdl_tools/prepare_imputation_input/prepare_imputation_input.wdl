task prepare_imputation_input{
    File bed_in
    File bim_in
    File fam_in
    File ref
    String? ref_group
    Float freq_diff_threshold = 0.2
    String? cohort

    String output_basename
    String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

    String docker = "rtibiocloud/prepare_imputation_input:v1.0_5a37578"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        set -e
        mkdir plink_input

        # Bed file preprocessing
        if [[ ${bed_in} =~ \.gz$ ]]; then
            # Append gz tag to let plink know its gzipped input
            gunzip -c ${bed_in} > plink_input/${input_prefix}.bed
        else
            # Otherwise just create softlink with normal
            ln -s ${bed_in} plink_input/${input_prefix}.bed
        fi

        # Bim file preprocessing
        if [[ ${bim_in} =~ \.gz$ ]]; then
            gunzip -c ${bim_in} > plink_input/${input_prefix}.bim
        else
            ln -s ${bim_in} plink_input/${input_prefix}.bim
        fi

        # Fam file preprocessing
        if [[ ${fam_in} =~ \.gz$ ]]; then
            gunzip -c ${fam_in} > plink_input/${input_prefix}.fam
        else
            ln -s ${fam_in} plink_input/${input_prefix}.fam
        fi

         python /opt/prepare_imputation_input.py \
            --bfile plink_input/${input_prefix} \
            --ref ${ref} \
            ${'--ref_group ' + ref_group} \
            --freq_diff_threshold ${freq_diff_threshold} \
            --out_prefix ${output_basename} \
            --working_dir $(pwd) \
            --plink /opt/plink \
            --bgzip /opt/bgzip \
            --bgzip_threads ${cpu} \
            ${'--cohort ' + cohort}
    >>>
    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        Array[File] vcf_out = glob("*.vcf.gz")
    }
}