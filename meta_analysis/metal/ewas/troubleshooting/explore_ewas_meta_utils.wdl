workflow wf { 
  Boolean comma_separated = true
  #File ewas_results = "Cannabis_GuLF_Model1_ANNOT_RESULTS_chr1_2021-06-09.csv"
  File ewas_results = "CCC.csv.gz"
  String study_basename = "gulf"

  Int probe_id_column = 1
  Int effect_size_column = 2
  Int standard_error_column = 3
  Int pvalue_column = 4
  Int chromosome_column = 6
  Array[Int] chromosomes_to_keep = [1]
  

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
  

  # Split file by chromosome
  call split_by_chromosome {
    input:
      comma_separated = comma_separated,
      ewas_results = to_split,
      chromosome_column = chromosome_column,
      study_basename = study_basename,
      chromosomes_to_keep = chromosomes_to_keep
  }

  Array[File] chr_split_results = split_by_chromosome.chr_files
  Array[Int] kept_chroms = split_by_chromosome.chr_order

  # keep only columns necessary for METAL
  scatter (chrom_order in range(length(chr_split_results))) {
    call keep_columns {
      input:
        infile =  chr_split_results[chrom_order],
        chromosome = kept_chroms[chrom_order],
        study_basename = split_by_chromosome.ewas_name,
        probe_id_column = probe_id_column,
        chromosome_column = chromosome_column,
        effect_size_column = effect_size_column,
        standard_error_column = standard_error_column,
        pvalue_column = pvalue_column
    }
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

  Boolean comma_separated
  File ewas_results 
  Int chromosome_column
  String study_basename 
  Array[Int] chromosomes_to_keep

  String docker = "ubuntu:18.04"
  Int cpu = 1
  Int mem = 2

  command <<<

  # format file
  # add option to take gzipped file
  # remove quotes and change to space separated

  if [[ ${comma_separated} == true ]]; then
    awk -F "," '
    {
      $1=$1
      gsub(/"/, "", $0)
    }1 
    ' OFS=" " ${ewas_results} > ${study_basename}_no_quote.txt

  else
    awk '
    {
      $1=$1
      gsub(/"/, "", $0)
    }1 
    ' OFS=" " ${ewas_results} > ${study_basename}_no_quote.txt

  fi
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

  #head ${study_basename}_no_quote_no_chr_prefix.txt

  # split results up by chromosome
  awk -v study=${study_basename} -v chr=${chromosome_column} \
    'NR==1{h = $0} NR>1{ print (!a[$chr]++ ? h ORS $0: $0) > study"_chr"$chr".txt"}' \
    ${study_basename}_no_quote_no_chr_prefix.txt

  # change chrX to chr23 if present
  mv ${study_basename}_chr[x,X].txt ${study_basename}_chr23.txt

  # keep only chromosomes specified
  mkdir keep_files/
  keep_chroms="${sep = ' ' chromosomes_to_keep}"

  for chr in $keep_chroms; do
    mv ${study_basename}_chr$chr.txt keep_files/
  done

  #rm ${study_basename}_chr*.txt # remove the chromosomes we don't want

  # get order of file:chromosome
  cd keep_files/
  ls ${study_basename}*t | perl -pe 's/.+chr//' \
    | perl -pe 's/.txt//'

  >>>

  output {
    Array[File] chr_files = glob("keep_files/${study_basename}*t")
    Array[String] chr_order = read_lines(stdout())
    String ewas_name = "${study_basename}"
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

task keep_columns {

  Int probe_id_column 
  Int chromosome_column
  Int effect_size_column
  Int standard_error_column
  Int pvalue_column

  Int chromosome
  String study_basename
  File infile

  String docker = "ubuntu:18.04"
  Int cpu = 1
  Int mem = 2
  
  command <<<

  # keep only specific columns
  echo -e "PROBE_ID CHR BETA SE P" > ${study_basename}_chr${chromosome}_specific_columns.tsv

    awk -v probe_id_column=${probe_id_column} \
        -v chromosome_column=${chromosome_column}  \
        -v effect_size_column=${effect_size_column} \
        -v pvalue_column=${pvalue_column} \
        -v standard_error_column=${standard_error_column} \
    '
        {print $probe_id_column, $chromosome_column, $effect_size_column, $standard_error_column, $pvalue_column}
    ' OFS=" " <(tail -n +2 ${infile}) >> ${study_basename}_chr${chromosome}_specific_columns.tsv

  >>>

  output {
  File metal_input = "${study_basename}_chr${chromosome}_specific_columns.tsv"
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }


  meta {
    description: "Keep only the specific columns needed for the metal analysis."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }

}

