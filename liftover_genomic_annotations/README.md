# LiftOver

Leverage the [`liftover`](https://neurogenomics.github.io/MungeSumstats/reference/liftover.html) function from the open-source package [MungeSumstats](https://github.com/neurogenomics/MungeSumstats) to seamlessly transfer genomic coordinates from one genome build to another.
This tool can convert GWAS summary statistics from GRCh38 to GRCh37 or vice versa.

We utilize their Docker image, available at [neurogenomicslab/mungesumstats](https://hub.docker.com/r/neurogenomicslab/mungesumstats), as a base image and include our custom Rscript that allows us to harness the `liftover` function within our WDL workflow.

## How to Run

Get started with the local runner and dev-toolkit [miniwdl](https://miniwdl.readthedocs.io/en/latest/index.html) in just a few steps:

1. Clone this repository.
1. Navigate to this folder.
1. Edit the `inputs.json` file to suit your needs.
1. Run the following command: `miniwdl zip main.wdl --input inputs.json`.
1. Then execute the workflow with `miniwdl run main.wdl.zip`.

## Inputs

**Description to be added.**


## Workflow Details

These WDL workflow files are specified to WDL syntax version 1.1.
We have configured the import statements to seamlessly integrate with both miniwdl and [AWS HealthOmics](https://aws.amazon.com/healthomics/).

## Contact

For any questions, suggestions, or further assistance, please feel free to reach out to the maintainer, Jesse Marks, at [jmarks@rti.org](mailto:jmarks@rti.org).
