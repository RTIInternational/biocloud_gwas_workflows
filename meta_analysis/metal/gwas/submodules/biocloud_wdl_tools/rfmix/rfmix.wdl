task rfmix{
    # Required stuff
    File query_vcf
    File ref_vcf
    File sample_map
    File genetic_map
    String output_basename
    String chr

    # Region to analyze <start_pos>-<end_pos>, in Mbp (decimal allowed)
    String? region

    # Random seed
    String? random_seed

    # Optional stuff
    Float? crf_spacing
    Float? rf_window_size
    Float? crf_weight
    Float? generations
    Int? em_iterations
    Boolean reanalyze_reference = false
    Int? node_size
    Int? trees
    Float? max_missing
    Int? bootstrap_mode
    Int? rf_min_snps
    String? debug_flag

    # Runtime environment
    String docker = "rtibiocloud/rfmix:v2_c39722e"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        set -e

        # Unzip query if needed (pigz is way faster)
        if [[ ${query_vcf} =~ \.gz$ ]]; then
            # Append gz tag to let plink know its gzipped input
            unpigz -p ${cpu} -c ${query_vcf} > query.vcf
        else
            # Otherwise just create softlink with normal
            ln -s ${query_vcf} query.vcf
        fi

        # Unzip ref if needed (pigz is way faster)
        if [[ ${ref_vcf} =~ \.gz$ ]]; then
            # Append gz tag to let plink know its gzipped input
            unpigz -p ${cpu} -c ${ref_vcf} > ref.vcf
        else
            # Otherwise just create softlink with normal
            ln -s ${ref_vcf} ref.vcf
        fi

        rfmix -f query.vcf \
            -r ref.vcf \
            -m ${sample_map} \
            -g ${genetic_map} \
            -o ${output_basename} \
            --chromosome=${chr} \
            ${'--crf-spacing=' + crf_spacing} \
            ${'--rf-window-size=' + rf_window_size} \
            ${'--crf-weight=' + crf_weight} \
            ${'--generations=' + generations} \
            ${'--em-iterations=' + em_iterations} \
            ${true='--reanalyze-reference' false='' reanalyze_reference} \
            ${'-n ' + node_size} \
            ${'-t ' + trees} \
            ${'--max-missing=' + max_missing} \
            ${'-b ' + bootstrap_mode} \
            ${'--rf-minimum-snps=' + rf_min_snps} \
            ${'--analyze-range=' + region} \
            ${'--debug=' + debug_flag} \
            ${'--n-threads=' + cpu} \
            ${'--random-seed=' + random_seed}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File msp_out = "${output_basename}.msp.tsv"
        File fb_out  = "${output_basename}.fb.tsv"
        File q_out   = "${output_basename}.rfmix.Q"
    }
}