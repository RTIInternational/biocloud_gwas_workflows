task terastructure_postprocess{
    File theta
    File fam
    File psam
    String ref_pop_type
    Array[String] ancestry_definitions
    String dataset_label
    String output_basename
    File? sample_study_groups

    # Runtime environment
    String docker = "rtibiocloud/terastructure_postprocessing:v1-37ef274"
    Int cpu = 2
    Int mem_gb = 4

    command {
        Rscript /opt/terastructure_postprocessing.R \
            --theta ${theta} \
            --fam ${fam} \
            --psam ${psam} \
            --ref_pop_type ${ref_pop_type} \
            --out_prefix ${output_basename} \
            ${'--sample_dataset_xref ' + sample_study_groups} \
            --ancestry_definition "${sep='" --ancestry_definition "' ancestry_definitions}"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[File] triangle_plots  = glob("${output_basename}*.png")
        Array[File] ancestry_thetas = glob("${output_basename}*theta.tsv")
        Array[File] ancestry_samples = glob("${output_basename}*keep")
        File misclassified_ref_samples = "${output_basename}_ref_misclassified.txt"
        File unclassified_samples = "${output_basename}_study_unclassified.txt"
    }
}

task order_by_ancestry{
    # Utility for aliging ancestry output files from terastructure_postprocess with a set of ancestry group names
    # This is kind of a hack because WDL doesn't yet support scattering over Map (hash) types
    # The fix here is to return matched lists of ancestry files and ancestry names so you can iterate over both
    # And know which ancestry file corresponds to a given ancestry
    Array[String] ancestry_files_in
    Array[String] ancestries

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        set -e

        # Loop through each ancesry and output filenames with the ancestry prefix
        for ancestry in ${sep=" " ancestries}
        do
            grep -i "$ancestry.keep" ${write_lines(ancestry_files_in)}

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

