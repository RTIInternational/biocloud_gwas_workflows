## GWAS Meta-analysis Workflow: METAL

This repository provides a workflow for running a Genome-Wide Association Study (GWAS) meta-analysis using the [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation) software tool. The workflow is written in the Workflow Description Language [(WDL)](https://github.com/openwdl/learn-wdl).


<br>

## Workflow Overview

This workflow utilizes Cromwell to execute the GWAS meta-analysis. While we provide a high-level overview here, detailed step-by-step instructions and video tutorials are available within the Omics Group team on Microsoft Teams:

* **Wiki:** Steps for running the workflow can be found on Teams: [HERE](https://researchtriangleinstitute.sharepoint.com/sites/OmicsGroup/_layouts/15/Doc.aspx?sourcedoc=%7Ba2b17bca-8f68-4450-a563-f80609bd497a%7D&action=edit&wd=target%28Computing%20Infrastructure.one%7C9f4912d8-44cd-4361-bebc-15e505bc676a%2FWDL%5C%2FCromwell%20Guide%7C80665feb-2dbf-481d-92d8-cf8c8e7d30dc%2F%29&wdorigin=703)  
* **Video Tutorial:**  A video tutorial for submitting WDL workflows to Cromwell can be found in the Files section of the Omics Group team on Microsoft Teams: [HERE](https://researchtriangleinstitute.sharepoint.com/sites/OmicsGroup/_layouts/15/stream.aspx?id=%2Fsites%2FOmicsGroup%2FShared%20Documents%2FComputing%20Infrastructure%2Fwdl%5Fsubmission%5Fwith%5Fcromwell%5Fserver%2Emp4&referrer=StreamWebApp%2EWeb&referrerScenario=AddressBarCopied%2Eview%2E9c6d32c0%2D893e%2D4391%2Db9db%2Df5658fc4896e)


User needs GWAS summary statitics with the following columns (headers can be named whatever):
   - rsID
   - chromosome
   - position
   - Beta (effect_size)
   - se_beta (standard error of beta coefficient)
   - A1 (coded allele)
   - A2 (noncoded allele)
   - P-value

For each cohort, the results should be combined into one input file — i.e., merge the individual chromosome results.

---

<br>

### `input.json` Configuration


The `input.json` file defines the parameters for the workflow. Here's a breakdown of the key parts:

- **container_source**: This specifies the Docker repository to pull images from (e.g., "ecr" or "dockerhub").
- **metal_gwas_meta_analysis_wf.study_basename**: This is a list of studies included in the meta-analysis. The order should match the order in subsequent entries.
- **metal_gwas_meta_analysis_wf.ancestry**: This defines the ancestry group for the meta-analysis (used for file naming).
- **metal_gwas_meta_analysis_wf.plot_basename**: This sets the base name for the generated QQ and Manhattan plots.
- **metal_gwas_meta_analysis_wf.full_results_name**: This sets the base name for the final results file.
- **metal_gwas_meta_analysis_wf.remove_singletons**: This defines whether to remove singletons (SNPs present in only one study). Set to "true" or "false".
- **metal_gwas_meta_analysis_wf.gwas_results**: This is a list of S3 paths to the GWAS results files for each study. Ensure the order matches the study list above.
- **metal_gwas_meta_analysis_wf.variant_id_column**: Column containing variant IDs (e.g., MarkerName) in each GWAS results file.
- **metal_gwas_meta_analysis_wf.chromosome_column**: Column containing chromosome information in each GWAS results file.  
- **metal_gwas_meta_analysis_wf.pos_column**: Column containing SNP position in each GWAS results file.
- **metal_gwas_meta_analysis_wf.coded_allele_column**: Column containing the coded allele in each GWAS results file.
- **metal_gwas_meta_analysis_wf.noncoded_allele_column**: Column containing the non-coded allele in each GWAS results file.
- **metal_gwas_meta_analysis_wf.effect_size_column**: Column containing the effect size in each GWAS results file.
- **metal_gwas_meta_analysis_wf.standard_error_column**: Column containing the standard error in each GWAS results file.
- **metal_gwas_meta_analysis_wf.pvalue_column**: Column containing the P-value in each GWAS results file.
- **metal_gwas_meta_analysis_wf.chromosomes_to_keep**: List of chromosomes to run the meta-analysis on. This allows you to run a meta based on a subset of chromosomes if desired.

\* The columns numbers are all 1-based instead of 0-based. So if the chromosome column is in the 2nd column, the column number is 2.
<br><br><br>

## To-do 
- [ ] Add 2df option
- [ ] Add option to apply [genomic control (GC)](https://en.wikipedia.org/wiki/Genomic_control).


# Authors
For any questions, please reach out to Jesse Marks (jmarks@rti.org).
