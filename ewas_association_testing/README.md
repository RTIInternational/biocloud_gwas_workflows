# EWAS
Details & instructions.

## Inputs

  "ewas_full.output_basename": "alspac_ea_ewas_model_1_results",
  "ewas_full.fdr_value": 0.10,
  "ewas_full.plot_colors": ["red", "blue"],
  "ewas_full.sample_name": "Sample_Name",
  "ewas_full.test_var": "cannabisUse",
  "ewas_full.covariates": ["age_at_DNAm", "sv1", "sv2", "sv3", "sv4", "sv5", "sv6", "Bcell", "CD4T", "CD8T", "Mono", "Gran", "NK"],
  "ewas_full.pheno_file": "s3://rti-cannabis/shared_data/post_qc/alspac/phenotypes/pheno_mothers_combined_FOM_TF1_3_n946_ewas_final.txt",
  "ewas_full.dnam_files": 

<details>
  <summary>docker</summary>

*String*<br>
  
Docker image containing the system dependencies and R packages necessary to run the workflow. Most up-to-date image as of 1-8-2021 is `ewas:v0.0.2_99db04b`. Visit rtibiocloud Docker Hub to get the latest version.
</details>
  
  
  
<details>
  <summary>fdr_value</summary>
  
  *Float*<br>
  
  False discovery rate.  
  </details>
  
  
  
  
<details>
  <summary></summary>
  
  
  </details>
  
  
  
  
  
<details>
  <summary></summary>
  
  
  </details>
  
  
  
  
  
  
<details>
  <summary></summary>
  
  
  </details>
  
  
  
  
  
  
<details>
  <summary></summary>
  
  
  </details>
  
  
  
  
  
<details>
  <summary></summary>
  
  
  </details>

  [ 

## Submitting the workflow
Microsoft Teams Omics Group [WDL/Cromwell Cheat Sheet](https://teams.microsoft.com/l/entity/com.microsoft.teamspace.tab.wiki/tab::61aecad5-13fa-4bde-adce-ba3b16950439?context=%7B%22subEntityId%22%3A%22%7B%5C%22pageId%5C%22%3A18%2C%5C%22origin%5C%22%3A2%7D%22%2C%22channelId%22%3A%2219%3Af42632e48b7c4b9e9f362afa1e4e1957%40thread.tacv2%22%7D&tenantId=2ffc2ede-4d44-4994-8082-487341fa43fb)
