task terastructure{
    File bed_in
    File bim_in
    File fam_in
    String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

    Int k
    Float rfreq_perc = 0.2
    String label = "terastruct"
    Boolean? force
    Int? seed
    Boolean? compute_beta
    File? id_map

    # Runtime environment
    String docker = "rtibiocloud/terastructure:v1.0-ee54b7a"
    Int cpu = 8
    Int mem_gb = 16

    command<<<
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

        # Count individuals and loci to be used
        num_indiv=$(wc -l plink_input/${input_prefix}.fam | cut -d" " -f 1)
        num_loci=$(wc -l plink_input/${input_prefix}.bim | cut -d" " -f 1)

        echo "Num indiv: $num_indiv"
        echo "Num loci: $num_loci"

        # Calculate report frequency as a percentage of the number of SNPs
        # This is basically just awk magic for CEIL(num_snps * rfreq_perc)
        rfreq=$(echo | awk -v var="$num_loci" '{print int(var*${rfreq_perc})}')

        echo "Rfreq: $rfreq"

        terastructure -file plink_input/${input_prefix}.bed \
            -n $num_indiv \
            -l $num_loci \
            -k ${k} \
            -rfreq $rfreq \
            ${'-label ' + label} \
            ${true='-force' false="" force} \
            ${'-idmap ' + id_map} \
            ${'-seed ' + seed} \
            ${true='-compute-beta' false="" compute_beta} \
            -nthreads ${cpu}

        # Move output files to working directory
        # Do this because terastruct makes random directory names with label appended to end
        mv ./*-${label}/theta.txt ./
        mv ./*-${label}/validation.txt ./
        mkdir logs
        mv ./*-${label}/ logs
        tar czvf terastructure_logs.tar.gz logs
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File admixture_proportions = "theta.txt"
        File validation = "validation.txt"
        File log_dir = "terastructure_logs.tar.gz"
    }

}