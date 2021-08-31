workflow wf {
  Boolean comma_separated = true
  #File ewas_results = "head_ewas_file.csv"
  Int chromosome_column = 6
  String study_basename = "gulf"
  File ewas_results = "Cannabis_GuLF_Model1_ANNOT_RESULTS_chr1_2021-06-09.csv"


  call split_by_chromosome {
    input:
      comma_separated = comma_separated,
      ewas_results = ewas_results,
      chromosome_column = chromosome_column,
      study_basename = study_basename
  }
}

task split_by_chromosome {

  Boolean comma_separated
  File ewas_results 
  Int chromosome_column
  String study_basename 
 # Array[Int] chromosomes_to_keep

  String docker = "ubuntu:18.04"
  Int cpu = 1
  Int mem = 2

  command <<<

  # format file
  # remove quotes and change to space separated
  if [[ ${comma_separated} == true ]]; then
    awk -F "," '
    {
      $1=$1
      gsub(/"/, "", $0)
    }1 
    ' OFS=" " ${ewas_results} > ${study_basename}_no_quote.txt

    # remove "chr" prefix if present from chromosome entries
    awk -v chr=${chromosome_column} '
        NR==1{print $0; next}
        NR>=2{
             chrom_prefix=substr($chr, 1, 3)
             chrom_suffix=substr($chr, 4)
             if (chrom_prefix == "chr" )
             {
                 $chr=chrom_suffix
             }
        }1
         ' ${study_basename}_no_quote.txt > ${study_basename}_no_quote_no_chr_prefix.txt
  else
    echo "Marks" 
  fi

  head ${study_basename}_no_quote_no_chr_prefix.txt


  >>>

  output {
   # String ewas_name = "${study_basename}"
   String out_string = stdout()
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }

  parameter_meta  {
    ewas_results: "The full path to the EWAS results going into the analysis."
    study_basename: "The name of the study/cohort going into the analysis. (e.g. uhs)"
    chromosome_column: "The 1-based chromosome column index. "
    chromosomes_to_keep: "Array of chromosomes that you want to analyze."
  }

  meta {
    description: "Input EWAS results with all chromosomes merged into one file. This task will then split these EWAS results up by chromosome. It will then remove the chromosomes not specified to keep. It also keeps track of which file is associated with which chromosome because the METAL results do not contain the chromosome information. Since we keep track of this, we can then go back and add in the chromosome column to the METAL results afterwards. We want the chromosome information because it is necessary for the Manhattan plot, among other reasons."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }
}



