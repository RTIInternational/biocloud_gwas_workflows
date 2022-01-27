task generate_gwas_plots{
    File summary_stats
    String col_id
    String col_chromosome
    String col_position
    String col_p
    String? variant_type_colname
    String output_basename

    Boolean? in_header = true
    Boolean? in_csv
    File? highlight_list

    Boolean? generate_manhattan_plot = true
    Boolean? generate_qq_plot = true
    Boolean? generate_snp_manhattan_plot
    Boolean? generate_indel_manhattan_plot
    Boolean? generate_snp_indel_manhattan_plot
    Boolean? generate_snp_qq_plot
    Boolean? generate_indel_qq_plot
    Boolean? generate_snp_indel_qq_plot

    Boolean? qq_lambda = true
    Boolean? qq_lines = true
    Boolean? qq_significance_line
    Int? manhattan_ylim
    Boolean? manhattan_no_line
    String? manhattan_odd_chr_color
    String? manhattan_even_chr_color
    String? manhattan_highlight_chr_color
    Float? manhattan_significance_value

    Array[String] cols_to_keep = select_all([col_id, col_chromosome, col_position, col_p, variant_type_colname])

    # Runtime options
    String docker = "rtibiocloud/generate_gwas_plots:v1_781a9d2"
    Int cpu = 1
    Int mem_gb = 16
    Int max_retries = 3

    command<<<
        set -e

        sumstats_file=${summary_stats}

        # unzip sumstats if gzipped so you can do some preprocessing
        if [[ "${summary_stats}" == *.gz ]]; then
            echo "un-gzipping sumstats..."
            gunzip -c ${summary_stats} > summary_stats.txt
            sumstats_file=summary_stats.txt
        fi

        # Get only the columns we actually need (reduce memory)
        echo 'BEGIN {
                split(cols,out,",")
              }
              NR==1 {
                for (i=1; i<=NF; i++)
                    ix[$i] = i
                for (i in out)
    	            printf "%s%s", $ix[out[i]], OFS
                print ""
              }
              NR>1 {
                for (i in out)
                    printf "%s%s", $ix[out[i]], OFS
                print ""
              }' > script.awk

        echo "Subsetting columns..."
        awk -f script.awk -v cols=${sep=',' cols_to_keep} $sumstats_file > fixed_sumstats.txt

        # Convert X chromosome to 23
        echo "Converting instances of X to 23..."
        sed -i 's/X/23/g' fixed_sumstats.txt

        /opt/generate_gwas_plots.R \
            --in fixed_sumstats.txt \
            --out ${output_basename} \
            --col_id ${col_id} \
            --col_chromosome ${col_chromosome} \
            --col_position ${col_position} \
            --col_p ${col_p} \
            ${'--col_variant_type ' + variant_type_colname} \
            ${true='--in_csv' false='' in_csv} \
            ${true='--in_header' false='' in_header} \
            ${'--highlight_list ' + highlight_list} \
            ${true='--generate_manhattan_plot' false='' generate_manhattan_plot} \
            ${true='--generate_qq_plot' false='' generate_qq_plot} \
            ${true='--generate_snp_manhattan_plot' false='' generate_snp_manhattan_plot} \
            ${true='--generate_snp_qq_plot' false='' generate_snp_qq_plot} \
            ${true='--generate_indel_manhattan_plot' false='' generate_indel_manhattan_plot} \
            ${true='--generate_indel_qq_plot' false='' generate_indel_qq_plot} \
            ${true='--generate_snp_indel_manhattan_plot' false='' generate_snp_indel_manhattan_plot} \
            ${true='--generate_snp_indel_qq_plot' false='' generate_snp_indel_qq_plot} \
            ${true='--qq_lambda' false='' qq_lambda} \
            ${true='--qq_lines' false='' qq_lines} \
            ${true='--qq_significance_line' false='' qq_significance_line} \
            ${true='--qq_significance_line' false='' qq_significance_line} \
            ${'--manhattan_ylim ' + manhattan_ylim} \
            ${'--manhattan_no_line ' + manhattan_no_line} \
            ${'--manhattan_odd_chr_color ' + manhattan_odd_chr_color} \
            ${'--manhattan_even_chr_color ' + manhattan_even_chr_color} \
            ${'--manhattan_highlight_chr_color ' + manhattan_highlight_chr_color} \
            ${'--manhattan_significance_value ' + manhattan_significance_value}
    >>>

    runtime{
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output {
        Array[File] plots = glob("${output_basename}.*.png")
   }

}
