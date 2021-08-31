workflow wf {
  File ewas_results = "CCC.csv.gz"

  # Unzip file if it needs to be unzipped
  if(basename(ewas_results) != basename(ewas_results,".gz")){
    call gunzip{
      input:
        in_file = ewas_results
    }
    File unzipped_ewas_results = gunzip.output_file
  }

  # This small section is because we have to choose whether to split the original input file or the unzipped input from previous task
  Array[File?] possible_files = [unzipped_ewas_results, ewas_results]
  File to_split = select_first(possible_files)

  call split_by_chromosome {
    input:
      ewas_results = to_split
  }
}



task gunzip {
  File in_file
  String out_filename = basename(in_file, ".gz")

  command{
      gunzip -c ${in_file} > ${out_filename}
  }
  output{
      File output_file = "${out_filename}"
  }
  runtime{
    docker: "ubuntu:18.04"
    cpu: "1"
    memory: "1 GB"
  }
}

task split_by_chromosome {

  File ewas_results 

  String docker = "ubuntu:18.04"
  Int cpu = 1
  Int mem = 2

  command <<<
  
  head ${ewas_results}

  >>>

  output {
    String outputt = stdout()
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }

  parameter_meta  {
    ewas_results: "The full path to the EWAS results going into the analysis."
    study_basename: "The name of the study/cohort going into the analysis. (e.g., uhs, gulf, or sister_study)"
    chromosome_column: "The 1-based chromosome column index."
    chromosomes_to_keep: "Array of chromosomes that you want to analyze."
  }

  meta {
    description: "Input EWAS results with all chromosomes merged into one file. This task will then split these EWAS results up by chromosome. It will then remove the chromosomes not specified to keep. It also keeps track of which file is associated with which chromosome because the METAL results do not contain the chromosome information. Since we keep track of this, we can then go back and add in the chromosome column to the METAL results afterwards. We want the chromosome information because it is necessary for the Manhattan plot, among other reasons."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }
}
#  if [[ ${study_basename: -2} == "gz" ]]; then
#    new_name=${study_basename:0:-4}
#  else
#    new_name=${study_basename}
#  fi
