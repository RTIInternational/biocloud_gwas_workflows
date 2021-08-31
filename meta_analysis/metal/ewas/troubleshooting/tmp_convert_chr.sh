#"cg25189904",-0.0254276827563594,0.00466681189772106,5.56553387048333e-08,2583,"chr1",68299493
#"cg12710152",0.000885355989119089,0.000193469092278691,4.9628709514271e-06,2583,"chr1",32716773
#"cg04259904",-0.00683236952625839,0.0015646473765048,1.31194036839774e-05,2583,"chr1",44399492
#"cg02818189",-0.010022727149719,0.00233907133315075,1.89629551284822e-05,2583,"chr1",25173236
#"cg13092108",-0.00607905917761719,0.00144572281510283,2.70310469089207e-05,2583,"chr1",26857284
#"cg10328583",-0.00392013768902231,0.000935335895581987,2.87047146437441e-05,2583,"chr1",6551073
#"cg03665853",0.00103727812988768,0.000248977043349762,3.20145376263081e-05,2583,"chr1",85219864
#"cg18992688",-0.00281072518706,0.000677703881043511,3.47301302230029e-05,2583,"chr1",206223241
#"cg13385220",-0.00925415984104166,0.00225290769811974,4.12376004183801e-05,2583,"chr1",202250453

cat head_ewas_file.csv

# remove quotes and change to space separated 
comma_separated="true"
if [[ ${comma_separated}=="true" ]]; then
    awk -F "," '
    {
    $1=$1
    gsub(/"/, "", $0) 
    }1   
    ' OFS=" " head_ewas_file.csv > head_ewas_file_no_quote.txt

    # remove "chr" prefix if present
    chr_col=6
    awk -v chr=${chr_col} '
        NR==1{print $0; next}
        NR>=2{
             chrom_prefix=substr($chr, 0, 3)
             chrom_suffix=substr($chr, 4)
             if ( chrom_prefix == "chr" )
             {
                 $chr=chrom_suffix
             }
        }1
         ' head_ewas_file_no_quote.txt # > head_ewas_file_no_quote_no_prefix.txt
fi
cat head_ewas_file_no_quote.txt


    
    chromosome_column=6
    # remove "chr" prefix if present from chromosome entries
    awk -v chr=${chromosome_column} '
        NR==1{print $0; next}
        NR>=2{
             chrom_prefix=substr($chr, 0, 3)
             chrom_suffix=substr($chr, 4)
             if (chrom_prefix == "chr" )
             {
                 $chr=chrom_suffix
             }
        }1
         ' head_ewas_file_no_quote.txt 



awk -v study=${study_basename} -v chr=${chromosome_column} \
  'NR==1{h = $0} NR>1{ print (!a[$chr]++ ? h ORS $0: $0) > study"_chr"$chr".txt"}' \
  <(zcat ${ewas_results})
