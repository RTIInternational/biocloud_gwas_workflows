# EWAS
Details & instructions.

## Inputs

<details>
  <summary> Click to Expand</summary>

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
  <summary>output_basename</summary>
  
  *String*<br>
  
A descriptive basename for the output file--e.g. "alspac_ea_ewas_model_1"
  </details>
  
  
  
  
  
<details>
  <summary>plot_colors</summary>
  
*Array[String]*<br>
  
An Array of two colors for the Manhattan plot. e.g. ["red", "blue"]
  </details>
  
  
  
  
<details>
  <summary>sample_name</summary>
 
  *String*<br>
  
 Sample name given in the DNAm data as well as the phenotype file.
</details>
  
  
  
  
  
  
<details>
  <summary>test_var</summary>
  
  *String*<br>
  
  Name of test variable.
  </details>
  
  
  
  
  
<details>
  <summary>covariates</summary>
  
  *Array[String]*
  
  List of covariates. 
  </details>
  

 
<details>
  <summary>pheno_file</summary>
  
  *String*<br>
  
  Location of phenotype file.
  </details>
  
  
  <details>
  <summary>dnam_files</summary>
  
  *Array[String]*<br>
  
  List of locations to the DNA methylation files.
  </details>
  
  </details>
  

## Submitting the workflow
Microsoft Teams Omics Group [WDL/Cromwell Cheat Sheet](https://teams.microsoft.com/l/entity/com.microsoft.teamspace.tab.wiki/tab::61aecad5-13fa-4bde-adce-ba3b16950439?context=%7B%22subEntityId%22%3A%22%7B%5C%22pageId%5C%22%3A18%2C%5C%22origin%5C%22%3A2%7D%22%2C%22channelId%22%3A%2219%3Af42632e48b7c4b9e9f362afa1e4e1957%40thread.tacv2%22%7D&tenantId=2ffc2ede-4d44-4994-8082-487341fa43fb)
