task terastructure_merger{
    File template_theta
    File template_fam
    File theta
    File fam
    String output_filename

    # Runtime environment
    String docker = "rtibiocloud/terastructure_merger:v1-77ee25f"
    Int cpu = 2
    Int mem_gb = 4

    command {
        Rscript /opt/terastructure_merger.R \
            --template_theta ${template_theta} \
            --template_fam ${template_fam} \
            --fam ${fam} \
            --theta ${theta} \
            --output_filename ${output_filename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File theta_out  = "${output_filename}"
    }
}