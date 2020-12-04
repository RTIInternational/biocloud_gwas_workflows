import "biocloud_gwas_workflows/association_testing/genesis/genesis_gwas_chr_wf.wdl" as GENESIS_CHR
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT

workflow genesis_gwas_wf{
    Array[File] genotype_files
    Array[File] output_file_prefixes
    Array[String] chrs
    File pheno_file
    String pheno_name
    # File? covar_file
    Array[String]? covars
    String geno_format
    String family
    String gxe

    Int genesis_cpu = 1
    Int genesis_mem_gb = 1

    # Split options
    Int? chunk_size
    String? variant_id_field
    Int? split_by_variant_cpu
    Int? split_by_variant_mem_gb

    # ID conversion parameters
    Array[File] id_ref_files
    String in_missing_allele
    String in_deletion_allele
    String ref_deletion_allele = "."
    String id_label = ""
    Int id_conversion_cpu = 1
    Int id_conversion_mem_gb = 1

    # TSV utils options
    Int tsv_append_mem_gb = 4

    # Do genesis chr workflow on each chromosome in parallel
    scatter(index in range(length(genotype_files))){

        call GENESIS_CHR.genesis_gwas_chr_wf as genesis{
            input:
                file_in_geno = genotype_files[index],
                file_out_prefix = output_file_prefixes[index],
                file_in_pheno = pheno_file,
                geno_format = geno_format,
                pheno = pheno_name,
                covars = covars,
                family = family,
                gxe = gxe,
                chr = chrs[index],
                genesis_cpu = genesis_cpu,
                genesis_mem_gb = genesis_mem_gb,
                chunk_size = chunk_size,
                variant_id_field = variant_id_field,
                split_by_variant_cpu = split_by_variant_cpu,
                split_by_variant_mem_gb = split_by_variant_mem_gb,
                tsv_append_mem_gb = tsv_append_mem_gb
        }

        # Convert ids in summary stats output to standard ids
        call IDCONVERT.convert_variant_ids{
            input:
                chr = chrs[index],
                in_file = genesis.summary_stats,
                in_header = 1,
                in_sep = "tab",
                ref = id_ref_files[index],
                in_id_col = 0,
                in_chr_col = 1,
                in_pos_col = 2,
                in_a1_col = 3,
                in_a2_col = 4,
                in_missing_allele = in_missing_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = ref_deletion_allele,
                output_filename = basename(genesis.summary_stats, ".tsv") + id_label + ".tsv.gz",
                output_compression = "gzip",
                cpu = id_conversion_cpu,
                mem_gb = id_conversion_mem_gb
        }

    }

    output{
        Array[File] unfiltered_sumstats = convert_variant_ids.output_file
    }
}