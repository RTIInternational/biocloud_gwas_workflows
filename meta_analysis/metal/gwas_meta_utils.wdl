task split_by_chromosome {

  String study_basename 
  File gwas_results 
  Array[Int] chromosomes_to_keep
  Int chromosome_column

  command <<<

  # split results up by chromosome
  awk -v study=${study_basename} -v chr=${chromosome_column} \
    'NR==1{h = $0} NR>1{ print (!a[$chr]++ ? h ORS $0: $0) > study"_chr"$chr".txt"}' \
    <(zcat ${gwas_results})
  
  # change chrX to chr23 if present
  mv ${study_basename}_chr[x,X].txt ${study_basename}_chr23.txt

  # keep only chromosomes specified
  mkdir keep_files
  keep_chroms="${sep = ' ' chromosomes_to_keep}"
  
  for chr in $keep_chroms; do
    mv ${study_basename}_chr$chr.txt keep_files
  done

  rm ${study_basename}_chr*.txt # remove the chromosomes we don't want

  # get order of file:chromosome
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
    docker: 'ubuntu:18.04'
  }

  parameter_meta  {
    gwas_results: "The full path to the GWAS results going into the analysis."
    study_basename: "The name of the study/cohort going into the analysis. (e.g. uhs)"
    chromosome_column: "The 1-based chromosome column index. "
    chromosomes_to_keep: "Array of chromosomes that you want to analyze."
  }

  meta {
    description: "Input GWAS results with all chromosomes merged into one file. This task will then split these GWAS results up by chromosome. It will then remove the chromosomes not specified to keep. It also keeps track of which file is associated with which chromosome because the METAL results do not contain the chromosome information. Since we keep track of this, we can then go back and add in the chromosome column to the METAL results afterwards. We want the chromosome information because it is necessary for the Manhattan plot, among other reasons."
    author: "Jesse Marks"
    email: "jmarks@rti.org"
  }
}


task keep_columns {

  Int variant_id_column 
  Int chromosome_column
  Int effect_size_column
  Int noncoded_allele_column
  Int coded_allele_column
  Int pvalue_column
  Int standard_error_column
  Int maf_column
  Int rsquared_column

  Int chromosome
  String study_basename
  File infile
  
  command <<<

  # keep only specific columns
  echo -e "VARIANT_ID\tCHR\tREF\tALT\tALT_EFFECT\tSE\tP\tMAF\tR_SQUARED" > head.txt

    cat head.txt > ${study_basename}_chr${chromosome}_specific_columns.tsv

    awk -v variant_id_column=${variant_id_column} \
        -v chromosome_column=${chromosome_column}  \
        -v noncoded_allele_column=${noncoded_allele_column} \
        -v coded_allele_column=${coded_allele_column} \
        -v effect_size_column=${effect_size_column} \
        -v pvalue_column=${pvalue_column} \
        -v standard_error_column=${standard_error_column} \
        -v maf_column=${maf_column} \
        -v rsquared_column=${rsquared_column} \
    '
        {print $variant_id_column, $chromosome_column, $noncoded_allele_column, $coded_allele_column, $effect_size_column, $standard_error_column, $pvalue_column, $maf_column, $rsquared_column}
    ' OFS="\t" <(tail -n +2 ${infile}) >> ${study_basename}_chr${chromosome}_specific_columns.tsv

  >>>


  output {
  File metal_input = "${study_basename}_chr${chromosome}_specific_columns.tsv"
  }
  
  runtime {
    docker: 'ubuntu:18.04'
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

  command <<<

    echo "
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
    docker: '404545384114.dkr.ecr.us-east-1.amazonaws.com/rtibiocloud/metal:v2020.05.05_1c7e830'
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
  
  command <<< 

    awk -v chrom=${chromosome} '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          f[$i] = i
        }
        print "Chromosome", "Position", $0
      }
      NR > 1 {
        {
          $2 = toupper($2);
          $3 = toupper($3);
          num_missing = gsub(/[?]/, "?", $f["Direction"]);
          num_cohorts = length($f["Direction"]);
          split($f["MarkerName"], pos, ":");
        }

        if ( num_cohorts == 1)
            print chrom, pos[2], $0;
        else if ( num_missing != (num_cohorts - 1) )
            print chrom, pos[2], $0;
      } ' OFS="\t" ${gwas_file} > chr${chromosome}_output_file.txt

  >>>

  output {
    File singletons_output = "chr${chromosome}_output_file.txt"
  }
  
  runtime {
    docker: 'ubuntu:18.04'
  }
  
  parameter_meta {
    chromosome: "Chromosome number."
    gwas_file: "Results file from METAL."
  }
  meta {
    description: "Remove singletons (SNPs that were only present in one study/cohort)"
    author: "Jesse Marks"
    email: "jmarks@rti.org"

  }

}

# make plot table
task merge_results {

  Array[File] gwas_results

  command <<<

    results_string="${sep=" " gwas_results}"

    python3 - << EOF
    import gzip
    import re

    results_list = "$results_string"
    outfile = "combined_sorted_results.txt"

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

    gzip combined_sorted_results.txt

  >>>
  
  output {
    File merged_results = "combined_sorted_results.txt.gz"
  }

  runtime {
    docker: "python:3.8" 
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
task final_results {
  File gwas_results
  Float pvalue

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
    File final_table = "final_results_p_lte_${pvalue}.tsv"
  }

  runtime {
    docker: "ubuntu:18.04" 
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
