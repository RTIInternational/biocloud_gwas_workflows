task cov_ldsc {
    File bed_in
    File bim_in
    File fam_in
    
    File cov_eigenvec
    String out_prefix

    String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

    # Runtime attributes
    String docker = "rtibiocloud/cov_ldsc:v1_78e7ecc"
    Int cpu
    Int mem_gb
    Int max_retries = 3

    command {
        
        # Get everything in the same directory
        mkdir plink_input

        # Bed file preprocessing
        ln -s ${bed_in} plink_input/${input_prefix}.bed

        # Bim file preprocessing
        ln -s ${bim_in} plink_input/${input_prefix}.bim

        # Fam file preprocessing
        ln -s ${fam_in} plink_input/${input_prefix}.fam

        # Run cov-LDSC python script
        set -e
        python /opt/ldsc.py \
            --bfile plink_input/${input_prefix} \
            --l2 \
            --ld-wind-cm 20 \
            --yes-really \
            --cov ${cov_eigenvec} \
            --out ${out_prefix} 
    }

    output {
        File m_File = "${out_prefix}.M"
        File m_5_50_File = "${out_prefix}.M_5_50"
        File ldscore_out = "${out_prefix}.ldscore.gz"
        File logFile = "${out_prefix}.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

}

