import "biocloud_gwas_workflows/biocloud_wdl_tools/genesis/genesis.wdl" as GENESIS
import "biocloud_gwas_workflows/biocloud_wdl_tools/gds_utils/gds_utils.wdl" as SPLIT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/helper_workflows/collect_large_file_list_wf.wdl" as COLLECT

workflow genesis_gwas_chr_wf{
    File file_in_geno
    File file_in_pheno
    File file_in_variant_list = "no_file"
    Boolean use_variant_list = false
    String file_out_prefix

    # Genesis options
    String geno_format       # Options: gds
    String pheno            # Column name in phenotype file
    Array[String]? covars   # Array of column names of covars
    String family           # Options: gaussian
    String? gxe        # Column name in phenotype file for GxE interaction
    String chr              # Chr being analyzed
    Int genesis_cpu = 1
    Int genesis_mem_gb = 1

    # Split options
    Boolean split_gds = false
    Int chunk_size = 500000
    String variant_id_field = "snp.id"
    Int split_by_variant_cpu = 1
    Int split_by_variant_mem_gb = 1

    # TSV utils options
    Int tsv_append_mem_gb = 4
    
    # Split genotypes for parallel processing
    if (split_gds) {
        call SPLIT.split_by_variant as split_by_variant{
            input:
                input_gds = file_in_geno,
                chunk_size = chunk_size,
                variant_id_field = variant_id_field,
                output_basename = "variant_list"
        }
    }

    # Loop through splits and do association testing on each
    if (defined(split_by_variant.split_lists) || use_variant_list) {
        Array[File] split_lists = select_first([split_by_variant.split_lists, [file_in_variant_list]])
        scatter(split_index in range(length(split_lists))){

            String split_file_out_prefix = if(split_gds) then (file_out_prefix + "_" + basename(split_lists[split_index], ".txt")) else file_out_prefix
            call GENESIS.genesis as split_genesis{
                input:
                    file_in_geno = file_in_geno,
                    file_in_pheno = file_in_pheno,
                    file_in_variant_list = split_lists[split_index],
                    file_out = split_file_out_prefix,
                    geno_format = geno_format,
                    pheno = pheno,
                    covars = covars,
                    family = family,
                    gxe = gxe,
                    chr = chr,
                    cpu = genesis_cpu,
                    mem_gb = genesis_mem_gb
            }
        }

        # Collect chunked sumstats files into single zip folder
        call COLLECT.collect_large_file_list_wf as collect_sumstats{
            input:
                input_files = split_genesis.sumstats_out,
                output_dir_name = file_out_prefix + "_genesis_output"
        }

        # Concat all sumstats files into single sumstat file
        call TSV.tsv_append as cat_sumstats{
            input:
                tsv_inputs_tarball = collect_sumstats.output_dir,
                output_filename = file_out_prefix + ".tsv",
                mem_gb = tsv_append_mem_gb
        }
    }

    if (!defined(split_by_variant.split_lists) && !use_variant_list) {
        call GENESIS.genesis{
            input:
                file_in_geno = file_in_geno,
                file_in_pheno = file_in_pheno,
                file_out = split_file_out_prefix,
                geno_format = geno_format,
                pheno = pheno,
                covars = covars,
                family = family,
                gxe = gxe,
                chr = chr,
                cpu = genesis_cpu,
                mem_gb = genesis_mem_gb
        }
    }

    File output_file = select_first([cat_sumstats.tsv_output, genesis.sumstats_out])

    output {
        File summary_stats = output_file
    }
}