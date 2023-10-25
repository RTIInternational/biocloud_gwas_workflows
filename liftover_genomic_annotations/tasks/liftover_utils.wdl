version 1.1

task liftOverTask {
  input {
    File file_name
    String output_name
    String sep
    String snp_name
    String chrom_name
    String pos_name
    String ref_genome
    String convert_ref_genome
    String chain_source

    String docker = "rtibiocloud/neurogenomicslab_mungesumstats:1.7.10"
    Int cpu = 1
    Int mem = cpu*2
  }

  command <<<
    Rscript /opt/neurogenomics_liftover.R \
        --file_name "~{file_name}" \
        --output_name "~{output_name}" \
        --sep "~{sep}" \
        --snp_name "~{snp_name}" \
        --chrom_name "~{chrom_name}" \
        --pos_name "~{pos_name}" \
        --ref_genome "~{ref_genome}" \
        --convert_ref_genome "~{convert_ref_genome}" \
        --chain_source "~{chain_source}"
  >>>

  output {
    File output_file = "~{output_name}"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{mem} GB"
  }
}
