import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK

workflow calculate_ld{

    Array[String] chrs

    # Input file parameters
    String input_format = "bed-bim-fam"
    Array[File] bed_files
    Array[File] bim_files
    Array[File] fam_files
    Array[File] vcf_files
    Array[File] bgen_files
    Array[File] sample_files

    # Output file parameters
    String output_basename
    String? output_format   # square, square0 triangle inter-chr
    Boolean? with_freqs
    Boolean? yes_really

    # Stat
    String? correlation_stat

    # LD window parameters
    Int? ld_window
    Int? ld_window_kb
    Float? ld_window_r2

    # Reference SNP parameters
    String? ld_snp
    File? ld_snp_list

    # D prime parameters
    String? dprime  # d, dprime, dprime-signed

    # Filtering options
    File? keep
    File? remove

    # Runtime/Resource options
    Int? cpu
    Int? mem_gb

    # Generate kinship matrix for each chromosome in parallel
    scatter(chr_index in range(length(chrs))){

        String chr = chrs[chr_index]

        call PLINK.calculate_ld as ld{
            input:
                input_format = input_format,
                bed = bed_files[chr_index],
                bim = bim_files[chr_index],
                fam = fam_files[chr_index],
                vcf = vcf_files[chr_index],
                bgen = bgen_files[chr_index],
                sample = sample_files[chr_index],
                output_basename = "${output_basename}_chr${chr}",
                output_format = output_format,
                with_freqs = with_freqs,
                yes_really = yes_really,
                correlation_stat = correlation_stat,
                ld_window = ld_window,
                ld_window_kb = ld_window_kb,
                ld_window_r2 = ld_window_r2,
                ld_snp = ld_snp,
                ld_snp_list = ld_snp_list,
                dprime = dprime,
                keep = keep,
                remove = remove,
                cpu = cpu,
                mem_gb = mem_gb
        }

    }

    output{
        Array[File] ld_files = ld.ld_file
        Array[File] log_files = ld.log_file
    }
}