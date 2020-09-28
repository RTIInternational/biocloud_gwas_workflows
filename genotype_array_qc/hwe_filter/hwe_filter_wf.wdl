import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/genotype_array_qc/normalize_sex_chr/normalize_sex_chr_wf.wdl" as NORM

workflow hwe_filter_wf{
    File bed_in
    File bim_in
    File fam_in
    String output_basename

    # HWE p-value filter threshold
    Float hwe_filter_pvalue
    String? hwe_mode

    # Whether or not to include non-founders (i.e. samples with parent_id set) in HWE calculations (default: False)
    Boolean? nonfounders

    # Workflow can optionally be used on a single chr
    Array[String]? chrs

    # Optional file of related individuals to exclude from HWE calculation
    File? related_samples

    Int plink_chr_cpu = 1
    Int plink_chr_mem_gb = 2
    Int plink_filter_cpu = 1
    Int plink_filter_mem_gb = 2

    # If user provides a custom chr list, use that instead to prevent having to merge-par unnecessarily
    # Default to merging sex chr (23 expected by norm_sex_wf)
    Array[String] expected_chrs = select_first([chrs, ["23"]])

    # Split PAR/NONPAR if not already split
    call NORM.normalize_sex_chr_wf{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            expected_chrs = ["23"],
            output_basename = "${output_basename}.sex_norm",
            build_code = "dummy",
            no_fail = true,
            plink_cpu = plink_filter_cpu,
            plink_mem_gb = plink_filter_mem_gb
    }

    File norm_bed = normalize_sex_chr_wf.bed_out
    File norm_bim = normalize_sex_chr_wf.bim_out
    File norm_fam = normalize_sex_chr_wf.fam_out

    # Get all chromosomes in bim if no list provided by user
    if(!defined(chrs)){
        call PLINK.get_bim_chrs{
            input:
                bim_in = norm_bim
        }
    }

    Array[String] scatter_chrs = select_first([chrs, get_bim_chrs.chrs])

    scatter(chr in scatter_chrs){

        # HWE filter on autosomes
        call PLINK.hardy as hwe{
            input:
                bed_in = norm_bed,
                bim_in = norm_bim,
                fam_in = norm_fam,
                chrs = [chr],
                filter_females = (chr == "23"),
                hwe_pvalue = hwe_filter_pvalue,
                hwe_mode = hwe_mode,
                nonfounders = nonfounders,
                remove_samples = related_samples,
                output_basename = "${output_basename}.chr.${chr}",
                cpu = plink_chr_cpu,
                mem_gb = plink_chr_mem_gb
        }
    }

    # Combine snps to remove across all chrs
    call UTILS.cat as cat_remove{
        input:
            input_files = hwe.remove,
            output_filename = "${output_basename}.hwe.remove"
    }

    # Remove failed HWE SNPs from full dataset
    call PLINK.make_bed as filter_hwe{
        input:
            bed_in = norm_bed,
            bim_in = norm_bim,
            fam_in = norm_fam,
            exclude = cat_remove.output_file,
            chrs = scatter_chrs,
            output_basename = "${output_basename}.hwe",
            cpu = plink_filter_cpu,
            mem_gb = plink_filter_mem_gb
    }

    output{
        File bed_out = filter_hwe.bed_out
        File bim_out = filter_hwe.bim_out
        File fam_out = filter_hwe.fam_out
    }

}