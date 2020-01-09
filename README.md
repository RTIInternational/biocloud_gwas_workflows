# biocloud_gwas_workflows
Portable WDL workflows for automating GWAS analysis pipelines at scale 

## Cloning this repo and biocloud_wdl_tools submodule together
	
	git clone --recurse-submodules https://github.rti.org/RTI/biocloud_gwas_workflows.git

## Pulling latest version of biocloud_wdl_tools into master branch

	git submodule update --remote biocloud_wdl_tools
 	git add biocloud_wdl_tools
	git commit -m "Pulled latest commit from biocloud_wdl_tools"
	git push origin master
