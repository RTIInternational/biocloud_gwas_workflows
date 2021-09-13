fileIn=/shared/jmarks/biocloud_gwas_workflows/meta_analysis/metal/ewas/cromwell-executions/metal_ewas_meta_analysis_wf/224bd07e-f7be-4e50-8fe4-84b1e65bf465/call-metal/shard-0/execution/afr_meta_analysis_chr1_output1.metal
# MarkerName      Allele1 Allele2 Effect  StdErr  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal
chromosome=1
#head  $fileIn

# format data
# add chromosome column

awk -v chrom=${chromosome} -F "\t" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
          f[$i] = i
        }
     print "Chromosome",$f["MarkerName"],$f["Effect"],$f["StdErr"],$f["P-value"],$f["Direction"]
    }

    NR > 1 {
        num_missing = gsub(/[?]/, "?", $f["Direction"])
        num_cohorts = length($f["Direction"])
        line = $f["MarkerName"]" "$f["Effect"]" "$f["StdErr"]" "$f["P-value"]" "$f["Direction"]
        #print num_cohorts 

        if(num_cohorts == 1) 
            print chrom, line
        else if( num_cohorts - num_missing > 1)
            print chrom, line
    }
'  $fileIn | head





awk -v chrom=${chromosome} -F "\t" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
              f[$i] = i
            }
        }
        { print "Chromosome",$f["MarkerName"],$f["Effect"],$f["StdErr"],$f["P-value"],$f["Direction"] }'

    {
        NR > 1 {
        {
          num_missing = gsub(/[?]/, "?", $f["Direction"])
          num_cohorts = length($f["Direction"])
          line = $f["MarkerName"]" "$f["Effect"]" "$f["StdErr"]" "$f["P-value"]" "$f["Direction"]
        }

        if ( num_cohorts == 1)
            {print chrom}
        else if ( num_cohorts - num_missing != 1)
            { print chrom}
        }
    } ' OFS=" " $fileIn | head
   
