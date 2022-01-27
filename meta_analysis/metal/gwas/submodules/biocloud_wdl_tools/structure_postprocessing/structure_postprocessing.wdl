task structure_postprocess{
    File parsed_structure_output
    File fam
    File psam
    String ref_pop_type
    Array[String] ancestry_definitions
    String output_basename

    # Optional values tat really
    File? sample_dataset_xref
    Boolean? color_by_dataset
    Float? tolerance

    # Runtime environment
    String docker = "rtibiocloud/structure_postprocessing:v1.0-05e9879"
    Int cpu = 1

    # This can sometimes get you into trouble but WTH we're gunna scale the memory based on the size of the files
    # It would be a shame if a really large job failed here because the user overlooked the amount of mem to this random R script
    Int mem_gb = ceil(size(parsed_structure_output, "GB") + size(fam, "GB") + size(psam, "GB")) + 1

    command {
        Rscript /opt/structure_postprocessing.R \
            --structure ${parsed_structure_output} \
            --fam ${fam} \
            --psam ${psam} \
            --ref_pop_type ${ref_pop_type} \
            --out_prefix ${output_basename} \
            ${'--tolerance ' + tolerance} \
            ${true="--color_by_dataset" false="" color_by_dataset} \
            ${'--sample_dataset_xref ' + sample_dataset_xref} \
            --ancestry_definition "${sep='" --ancestry_definition "' ancestry_definitions}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[File] triangle_plots  = glob("${output_basename}*.png")
        Array[File] ancestry_proportions = glob("${output_basename}*admix.tsv")
        Array[File] ancestry_samples = glob("${output_basename}*keep")
        File unclassified_samples = "${output_basename}_unclassified.txt"
    }
}

task order_by_ancestry{
    # Utility for aliging ancestry output files from terastructure_postprocess with a set of ancestry group names
    # This is kind of a hack because WDL doesn't yet support scattering over Map (hash) types
    # The fix here is to return matched lists of ancestry files and ancestry names so you can iterate over both
    # And know which ancestry file corresponds to a given ancestry
    Array[String] ancestry_files_in
    Array[String] ancestries
    String? input_basename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        set -e

        # Loop through each ancesry and output filenames with the ancestry prefix
        for ancestry in ${sep=" " ancestries}
        do
            grep -i "${input_basename}_$ancestry.keep" ${write_lines(ancestry_files_in)}

        done
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[String] ancestry_file_out = read_lines(stdout())
    }
}

