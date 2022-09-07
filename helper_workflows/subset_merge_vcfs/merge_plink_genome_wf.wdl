import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

task merge_plink{
    # This merging is specific to TOPMed whole genome sequencing (WGS) data.

    Array[File] bed_in
    Array[File] bim_in
    Array[File] fam_in
    String output_basename
    
    String input_prefix = basename(sub(fam_in, "\\.gz$". ""), ".fam")

    String docker = "rtibiocloud/plink:v1.9_178bb91"
    Int cpu
    Int mem_gb
    Int max_retries = 3

    command <<
        mkdir plink_input

        # Bed file preprocessing
        for file in ${sep=" " bed_in}; do
            String input_prefix_bed = basename(sub($file, "\\.gz$". ""), ".bed")
            if [[ ${file} =~ \.gz$ ]]; then
                # Append gz tag to let plink know its gzipped input
                unpigz -p ${cpu} -c ${file} > plink_input/${input_prefix}.bed
            else
                # Otherwise just create softlink with normal
                ln -s ${file} plink_input/${input_prefix}.bed
            fi
            echo "plink_input/${input_prefix}.bed" >> plink_input/bed_files.txt
        done

        # Bim file preprocessing
        for file2 in ${sep=" " bim_in}; do
            if [[ ${file2} =~ \.gz$ ]]; then
                unpigz -p ${cpu} -c ${file2} > plink_input/${input_prefix}.bim
            else
                ln -s ${file2} plink_input/${input_prefix}.bim
            fi
            echo "plink_input/${input_prefix}.bim" >> plink_input/bim_files.txt
        done

        # Fam file preprocessing
        for file3 in ${sep=" " fam_in}; do
            if [[ ${fam_in} =~ \.gz$ ]]; then
                unpigz -p ${cpu} -c ${fam_in} > plink_input/${input_prefix}.fam
            else
                ln -s ${fam_in} plink_input/${input_prefix}.fam
            fi
            echo "plink_input/${input_prefix}.fam" >> plink_input/fam_files.txt
        done

        # Merge bed/bim/fam links into merge-list file
        paste -d " " plink_input/bed_files.txt plink_input/bim_files.txt plink_input/fam_files.txt > plink_input/merge_list.txt

        # Now run plink
        plink --make-bed \
            --threads ${cpu} \
            --set-missing-var-ids @:#:\$r:$a \
            --new-id-max-allele-len 300 \
            --max-alleles 2 \
            --merge-list plink_input/merge_list.txt \
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

workflow merge_plink_wf{
    Array[File] beds_in
    Array[File] bims_in
    Array[File] fams_in

    String output_basename
    Int cpu
    Int mem_gb

    call merge_plink{
        input:
            bed_in = beds_in
            bim_in = bims_in,
            fam_in = fams_in,
            output_basename = output_basename

            cpu = cpu
            mem_gb = mem_gb
    }
    
    output {
        File bed = merge_plink.bed_out
        File bim = merge_plink.bim_out
        File fam = merge_plink.fam_out
        File log = merge_plink.plink_log
    }
    
}

