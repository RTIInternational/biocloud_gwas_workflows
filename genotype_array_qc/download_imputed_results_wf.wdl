task download_imputed_vcf{
    String url
    String pwd
    String output_prefix

    # Runtime environment
    String docker = "rtibiocloud/plink:v1.9_178bb91"
    Int cpu = 2
    Int mem_gb = 4

    command <<<
        set -e

        # Download from imputation server
        wget ${url}

        # Unzip to get vcf and info files
        unzip -P ${pwd} $(ls *.zip)

        # Add output prefix to files
        mv $(ls *.vcf.gz) "${output_prefix}.$(ls *.vcf.gz)"
        mv $(ls *.info.gz) "${output_prefix}.$(ls *.info.gz)"
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[File] vcf_out = glob("${output_prefix}*.vcf.gz")
        Array[File] info_out = glob("${output_prefix}*.info.gz")
    }
}

workflow download_imputed_results_wf{
    Array[String] urls
    String pwd
    String output_prefix

    # Scatter across cohorts
    scatter(url in urls){
        call download_imputed_vcf{
            input:
                url = url,
                pwd = pwd,
                output_prefix = output_prefix
        }
    }

    output{
        Array[Array[File]] vcf = download_imputed_vcf.vcf_out
        Array[Array[File]] info = download_imputed_vcf.info_out
    }
}