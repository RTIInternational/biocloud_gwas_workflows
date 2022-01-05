# LiftOver
LiftOver is performed to bring all your genetical analyses to the same reference build.
Use this workflow to perform a genome-wide liftover on your GWAS summary statitics.

For example, convert GWAS summary statistics from GRCh38 to GRCh37, or vice versa.

# Step 1: Select chainfile
* Navigate to the Human genome liftOver files at http://hgdownload.soe.ucsc.edu/downloads.html#human.
* Click "LiftOver files" under the sub-header of your starting genome build. For example, if my summary statistics are currently in build 38 I would click on [LiftOver files](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) under the sub-header `Dec. 2013 (GRCh38/hg38)`.
* Copy the link to the desired build conversion. For example, if my sumstats are currenlty in build 38 and I want to convert to build 37 the link would be `http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz`.

# Step 2: Prepare input file
The inputs JSON file is where the parameters are set. 
<details>
  <summary>genome_liftover.input_sumstats</summary>
  
  S3 Path to GWAS summary statistics.
  * The results must be merged into one input file. If your results are split by chromosome, you will need to merge them before performing a genome-wide liftover. Though, the workflow will still work if you only want to perform the liftover on one set of chromosome results. 
  * Make sure the results in the S3 Standard storage tier and not archived. They must be restored if they are 
  </details>
  
  
  
  
  <details>
  <summary>genome_liftover.final_file</summary>
  
  The desired name of output file.
  </details>
  
  
  
  
  
  <details>
  <summary>genome_liftover.chainfile</summary>
  
  Paste the link to the chain file you copied in Step 1.
  </details>
  
  
  
  
  <details>
  <summary>genome_liftover.chromosome_col</summary>
  
   Zero-based
  </details>
  
  
  
  <details>
  <summary>genome_liftover.position_col</summary>
  
  Zero-based
  </details>
