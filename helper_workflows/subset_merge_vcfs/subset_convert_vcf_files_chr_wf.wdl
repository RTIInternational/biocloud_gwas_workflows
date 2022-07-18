import "biocloud_gwas_workflows/biocloud_wdl_tools/bcftools/bcftools.wdl" as SUBSET

task convert_vcf_to_bed{
    # This conversion is specific to TOPMed whole genome sequencing (WGS) data. It includes filtering for (1) only biallelic SNPs, and (2) only SNPs which pass the TOPMed SVM filter.

    File vcf_in
    String output_basename
    String input_prefix = basename(sub(vcf_in, "\\.gz$", ""), ".vcf")

    String docker = "rtibiocloud/plink:v2.0_888cf13"
    Int cpu = 2
    Int mem_gb = 2
    Int max_retries = 3

    command <<<
        set -e
        mkdir plink_input

        # VCF file preprocessing
        if [[ ${vcf_in} =~ \.gz$ ]]; then
            # Unzip
            unpigz -p ${cpu} -c ${vcf_in} > plink_input/${input_prefix}.vcf
        else
            # Otherwise just create softlink with normal
            ln -s ${vcf_in} plink_input/${input_prefix}.vcf

        # Convert
        plink2 \
            --vcf plink_input/${input_prefix}.vcf \
            --vcf-idspace-to_ --const-fid --allow-extra-chr 0 \
            --set-missing-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 23 missing \
            --vcf-filter PASS \
            --max-alleles 2 \
            --make-bed \
            --out ${output_basename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File bed_out = "${output_basename}.bed"
        File bim_out = "${output_basename}.bim"
        File fam_out = "${output_basename}.fam"
        File plink_log = "${output_basename}.log"
    }
}

workflow subset_convert_vcf_chr_wf{
    # BCFtools view options
    File vcf_in
    File samples_file
    String output_filename
    String output_type

    String output_basename

    # Do subset workflow
    call SUBSET.view as view{
        input:
            vcf_in = vcf_in,
            samples_file = samples_file,
            output_filename = output_filename,
            output_type = output_type
    }

    call convert_vcf_to_plink{
        input:
            vcf_in = view.vcf_out,
            output_basename = output_basename
    }
    
    output {
        File convert.bed = convert.bed_out
        File convert.bim = convert.bim_out
        File convert.fam = convert.fam_out
        File convert.log = convert.plink_log
    }
}

