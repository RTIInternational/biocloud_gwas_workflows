task gunzip {
  File in_file 
  # this should probs be in {} since it internal var
  String out_filename = basename(in_file, ".gz")

  String docker = "ubuntu:22.04"
  Int cpu = 1
  Int mem = 2


  command{
      gunzip -c ${in_file} > ${out_filename}
  }
  output{
      File output_file = "${out_filename}"
  }
    runtime{
      docker: docker
      cpu: cpu
      memory: "${mem} GB"
    }
}


task split_by_chromosome {

  Boolean comma_separated
  File gwas_results 
  Int chromosome_column
  String study_basename 
  Array[Int] chromosomes_to_keep

  String docker = "ubuntu:22.04"
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
    ' OFS=" " ${gwas_results} > ${study_basename}_no_quote.txt
  else
    awk '
    {
      $1=$1
      gsub(/"/, "", $0)
    }1 
    ' OFS=" " ${gwas_results} > ${study_basename}_no_quote.txt

  fi

  # remove "chr" prefix if present in the chromosome entries
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

  # get the order of file:chromosome
  cd keep_files/
  ls ${study_basename}*t | perl -pe 's/.+chr//' \
    | perl -pe 's/.txt//'

  >>>

  output {
    Array[File] chr_files = glob("keep_files/${study_basename}*t")
    Array[String] chr_order = read_lines(stdout())
    String gwas_name = "${study_basename}"
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }

  parameter_meta  {
    gwas_results: "The full path to the GWAS results going into the analysis."
    study_basename: "The name of the study/cohort going into the analysis. (e.g., uhs, gulf, or sister_study)"
    chromosome_column: "The 1-based chromosome column index."
    chromosomes_to_keep: "Array of chromosomes that you want to analyze."
  }

  meta {
    description: "Input GWAS results with all chromosomes merged into one file. This task will then split these GWAS results up by chromosome. The first step is to remove any quotes from the file and change to space separated. It will then remove the chromosomes not specified to keep. It also keeps track of which file is associated with which chromosome because the METAL results do not contain the chromosome information. Since we keep track of this, we can then go back and add in the chromosome column to the METAL results afterwards. We want the chromosome information because it is necessary for the Manhattan plot, among other reasons."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }
}


task keep_columns {

  Int variant_id_column 
  Int chromosome_column
  Int pos_column
  Int noncoded_allele_column
  Int coded_allele_column
  Int effect_size_column
  Int standard_error_column
  Int pvalue_column

  Int chromosome
  String study_basename
  File infile

  String docker = "ubuntu:22.04"
  Int cpu = 1
  Int mem = 2
  
  command <<<

  # keep only specific columns
  echo -e "VARIANT_ID CHROMOSOME POSITION REF ALT ALT_EFFECT SE P" > ${study_basename}_chr${chromosome}_specific_columns.txt

    awk -v variant_id_column=${variant_id_column} \
        -v chromosome_column=${chromosome_column}  \
        -v pos_column=${pos_column} \
        -v noncoded_allele_column=${noncoded_allele_column} \
        -v coded_allele_column=${coded_allele_column} \
        -v effect_size_column=${effect_size_column} \
        -v standard_error_column=${standard_error_column} \
        -v pvalue_column=${pvalue_column} \
    '
        {print $variant_id_column, $chromosome_column, $pos_column, $noncoded_allele_column, $coded_allele_column, $effect_size_column, $standard_error_column, $pvalue_column}
    ' OFS=" " <(tail -n +2 ${infile}) >> ${study_basename}_chr${chromosome}_specific_columns.txt

  >>>

  output {
  File metal_input = "${study_basename}_chr${chromosome}_specific_columns.txt"
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

task run_metal {
  String ancestry 
  Int chromosome
  Array[File] gwas_files
  Int cpu = 1
  Int mem = 4
  String docker

  command <<<

    echo "
    TRACKPOSITIONS ON
    SCHEME STDERR
    PVALUE P
    MARKER VARIANT_ID
    ALLELE ALT REF
    EFFECT ALT_EFFECT
    STDERR SE
    GENOMICCONTROL OFF 
    PROCESS ${sep= ' \n PROCESS ' gwas_files}
    OUTFILE ${ancestry}_meta_analysis_chr${chromosome}_output .metal
    ANALYZE HETEROGENEITY
    QUIT" > metal.txt

    # execute metal file
    /opt/metal metal.txt

  >>>

  output {
    File metal_results = "${ancestry}_meta_analysis_chr${chromosome}_output1.metal"
 }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }

  parameter_meta  {
    ancestry: "Ancestry of study/cohort (e.g. afr, amr, or eur)"
    chromosome: "Chromosome number of input file."
    gwas_files: "Array of input GWAS summary stats for the specified chromosome."
  }

  meta {
    description: "Perform meta-analysis on the input files of the specified chromosome."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }
}


task exclude_singletons {
  Int chromosome
  File gwas_file
  String remove_singletons
  String docker = "ubuntu:22.04"
  Int cpu = 1
  Int mem = 2
  
  command <<< 

    # format data
      awk -v decision="${remove_singletons}" '
    # print the header
    NR == 1 {
      # redundant asignment so we can asign the OFS
      $1 = $1 
        OFS=" "
      for (i = 1; i <= NF; i++) {
        f[$i] = i
      }
      print $0
    }

    # print pruned line if it is not a singleton
    NR > 1 {
        if (decision == "true") {
          $f["Allele1"] = toupper($f["Allele1"])
          $f["Allele2"] = toupper($f["Allele2"])
          num_missing = gsub(/[?]/, "?", $f["Direction"])
          num_cohorts = length($f["Direction"])
          #line = $f["Chromosome"]" "$f["Position"]" "$f["MarkerName"]" $f["Allele1"] $f["Allele2"] "$f["Effect"]" "$f["StdErr"]" "$f["P-value"]" "$f["Direction"]

          if(num_cohorts == 1)
              print $0
              #print line
          else if( num_cohorts - num_missing > 1 )
              print $0
              #print line
        }
        
        else {
          $f["Allele1"] = toupper($f["Allele1"])
          $f["Allele2"] = toupper($f["Allele2"])
          print $0
        }
    } ' ${gwas_file} > chr${chromosome}_output_file.txt
  >>>

  output {
    File singletons_output = "chr${chromosome}_output_file.txt"
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }
  
  parameter_meta {
    gwas_file: "Results file from METAL."
    chromosome: "The chromosome number (Int)."
    remove_singletons: "String <true> or <false> that determines if we will exclude singletons."
  }
  meta {
    description: "Capitalize alleles from METAL output and, if specified, remove singletons (variants/SNPs that were only present in one study/cohort)."
    author: "Jesse Marks"
    email: "jmarks@rti.org"

  }
}

# make plot table
task merge_final_results {

  Array[File] gwas_results
  String full_results_name = "final_results_table"

  String docker = "python:3.12-alpine"
  Int cpu = 1
  Int mem = 2

  command <<<

    results_string="${sep=" " gwas_results}"

    python3 - << EOF
    import gzip
    import re

    results_list = "$results_string"
    outfile = "${full_results_name}.txt"

    results_list = results_list.split()

    chroms_dict = {}
    # populate dictionary with chrom (key) and path to results (value)
    for results in results_list:
        substring = re.sub(r'.+\/chr', '', results)
        chrom = re.sub(r'_output_file.txt', '', substring)
        chroms_dict[int(chrom)] = results

    # combine sorted chromosome results files
    with open(outfile, 'w') as outF:
        count = 0 # print the header line if it's the first file

        for chrom in sorted(chroms_dict):
            infile = str(chroms_dict[chrom])

            with open(infile) as inF:
                head = inF.readline()
                if count == 0:
                    outF.write(head)

                # append the rest of the lines
                line = inF.readline()
                while line:
                    outF.write(line)
                    line = inF.readline()

            count += 1

    EOF

    gzip ${full_results_name}.txt

  >>>
  
  output {
    File merged_results = "${full_results_name}.txt.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }

  parameter_meta {
    gwas_results: "Array of results split by chromosome."
  }
  
  meta {
    description: "Merge the chromosome results back together so they can be input into the gwas plotting workflow."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }
}

# make pvalue filtered file we will put default at Pvalue<0.001 (user can specify though)
# go ahead and order this too
task top_results {
  File gwas_results
  Float pvalue

  String docker = "ubuntu:22.04"
  Int cpu = 1
  Int mem = 2

  command <<<

    awk -v pval="${pvalue}" '
      NR==1 {
        for (i=1; i<=NF; i++)
          F[$i] = i;
          print $0;
      }
      NR>1 {
        if ($F["P-value"] <= pval)
          print $0;
      }' <(zcat ${gwas_results}) > "final_results_p_lte_${pvalue}.tsv"

  >>>

  output {
    File final_results_pvalue_filtered = "final_results_p_lte_${pvalue}.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "${mem} GB"
  }

  parameter_meta {
    gwas_results: "The GWAS meta-analysis results containing only the chromosomes of interest."  
    pvalue: "P-value threshold to filter by."
  }

  meta {
    description: "Create the final, P-value filtered file."
    author: "Jesse Marks"
    email: "jmarks@rti.org"

  }
}        
