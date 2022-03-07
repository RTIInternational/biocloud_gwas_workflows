import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_variant_ids/convert_variant_ids.wdl" as IDCONVERT
import "biocloud_gwas_workflows/biocloud_wdl_tools/tsv_utils/tsv_utils.wdl" as TSV
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS
import "biocloud_gwas_workflows/genotype_array_qc/normalize_sex_chr/normalize_sex_chr_wf.wdl" as NORM

task chr_sort_check{
    # Return whether input chrs are in numerical sort order
    Array[String] chrs
    String docker = "ubuntu:18.04"

    command<<<
        sort -n ${write_lines(chrs)} > sorted.txt
        diff sorted.txt ${write_lines(chrs)} > diff.txt
        if [ $? -eq 0 ];then
            echo "false"
        else
            echo "true"
        fi
    >>>

    runtime {
        docker: docker
        cpu: 1
        memory: "100 MB"
    }

    output{
        Boolean failed = read_boolean(stdout())
    }

}

task label_duplicate_variants{
    File bim_in
    String output_basename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command <<<
        set -e
        # Get sorted list of all variants
        perl -lane 'print $F[1];' ${bim_in} | sort > tmp.variants.sorted

        # Get sorted list of unique variants
        sort -u tmp.variants.sorted > tmp.variants.sorted.unique

        # Get list of duplicate variants
        comm -23 tmp.variants.sorted \
            tmp.variants.sorted.unique | sort -u > tmp.variants.duplicates

        # Append ___X (where X is a unique number for the variant) to the end of the variant IDs for duplicates
        perl -lane '
            BEGIN {
                %duplicates = ();
                open(DUPLICATES, "tmp.variants.duplicates");
                while (<DUPLICATES>) {
                    chomp;
                    $duplicates{$_} = 1;
                }
                close DUPLICATES;
            }
            if (exists($duplicates{$F[1]})) {
                $F[1] = $F[1]."___".($duplicates{$F[1]}++);
            }
            print join("\t", @F);' ${bim_in} > ${output_basename}.bim
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File bim_out = "${output_basename}.bim"
    }
}

task get_duplicate_variant_ids{
    File bim_in
    String output_filename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        grep ___ ${bim_in} > ${output_filename}
        touch ${output_filename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File duplicate_ids = "${output_filename}"
    }
}

task get_variants_to_remove{
    File bed_in
    File bim_in
    File fam_in
    File duplicate_snp_ids
    String? chr
    String output_basename
    String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

    # Runtime environment
    String docker = "rtibiocloud/plink:v1.9_178bb91"
    Int cpu = 4
    Int mem_gb = 8

    command <<<
        set -e

        mkdir plink_input

        # Bed file preprocessing
        if [[ ${bed_in} =~ \.gz$ ]]; then
            # Append gz tag to let plink know its gzipped input
            unpigz -p ${cpu} -c ${bed_in} > plink_input/${input_prefix}.bed
        else
            # Otherwise just create softlink with normal
            ln -s ${bed_in} plink_input/${input_prefix}.bed
        fi

        # Bim file preprocessing
        if [[ ${bim_in} =~ \.gz$ ]]; then
            unpigz -p ${cpu} -c ${bim_in} > plink_input/${input_prefix}.bim
        else
            ln -s ${bim_in} plink_input/${input_prefix}.bim
        fi

        # Fam file preprocessing
        if [[ ${fam_in} =~ \.gz$ ]]; then
            unpigz -p ${cpu} -c ${fam_in} > plink_input/${input_prefix}.fam
        else
            ln -s ${fam_in} plink_input/${input_prefix}.fam
        fi

        # Get call rates of duplicate SNPs
        plink --bfile plink_input/${input_prefix} \
            --missing \
            --extract ${duplicate_snp_ids} \
            ${'--chr ' + chr} \
            --out ${output_basename}

        # Create remove list containing dupliicate with higher missing rate
        tail -n +2 ${output_basename}.lmiss |
        perl -lane '
            BEGIN {
                %missingness = ();
                %extract = ();
                @variants = ();
            }
            push(@variants, $F[1]);
            $F[1] =~ /^(\S+)___\d+/;
            $baseName = $1;
            if (exists($missingness{$baseName})) {
                if ($F[4] < $missingness{$baseName}) {
                $missingness{$baseName} = $F[4];
                $extract{$baseName} = $F[1];
                }
            } else {
                $missingness{$baseName} = $F[4];
                $extract{$baseName} = $F[1];
            }
            END {
                foreach $variant (@variants) {
                    $variant =~ /^(\S+)___\d+/;
                    $baseName = $1;
                    if ($extract{$baseName} ne $variant) {
                        print $variant;
                    }
                }
            }' > ${output_basename}.duplicates_to_remove.txt
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File dup_missing_report = "${output_basename}.lmiss"
        File dups_to_remove = "${output_basename}.duplicates_to_remove.txt"
    }
}

task fix_ids{
    File bim_in
    String output_basename

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        set -e
        perl -lne 's/___\d+//; print;' ${bim_in} > ${output_basename}.bim
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File bim_out = "${output_basename}.bim"
    }
}

workflow impute2_id_conversion_wf{
    File bed_in
    File bim_in
    File fam_in
    Array[String] chrs
    Array[File] id_legend_files
    Boolean remove_duplicates = true
    String output_basename

    # PAR/NONPAR Split/Merge parameters
    String build_code
    Boolean no_fail = true
    String in_monomorphic_allele
    String in_deletion_allele
    String ref_deletion_allele = "."
    Boolean rescue_rsids = true

    # Resources subtasks
    # The only one you'll likely need to play for huge files is the id_convert cpu/mem
    # But mainly you just need to make sure the mem is roughly the size of the largest id_legend file
    Int plink_cpu = 1
    Int plink_mem_gb = 2

    Int id_convert_cpu = 1
    Int id_convert_mem_gb = 3

    Int duplicate_id_cpu = 2
    Int duplicate_id_mem_gb = 6
    
    String docker_ubuntu = "404545384114.dkr.ecr.us-east-1.amazonaws.com/ubuntu:18.04"

    # Make sure chromsomes are provided in numerical sort order and error out if they aren't
    # This is to ensure that the id_legend_files, chrs, and bim file all have same sort order
    # Necessary because we can then operate directly on the bim file without having to make chr-level that need to be merged at the end (costly, slow)
    call chr_sort_check as check_id_legend_chr_sort{
        input:
            chrs = chrs,
            docker = docker_ubuntu
    }

    if(check_id_legend_chr_sort.failed){
        call UTILS.raise_error{
            input:
                msg = "Input chromosomes MUST be provided in numerical sort order!",
                docker = docker_ubuntu
        }
    }

    # Conditionally split/merge PAR/NONPAR regions of input chrX based on list of provided chromosomes
    call NORM.normalize_sex_chr_wf{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            expected_chrs = chrs,
            build_code = build_code,
            no_fail = no_fail,
            output_basename = output_basename,
            plink_cpu = plink_cpu,
            plink_mem_gb = plink_mem_gb
    }

    # Need to subset to make sure only chrs are included in bed file so we can work with the bim by itself
    # Have to use Plink1.9 here bc it sorts automatically and plink2 doesn't sort at all (annoying)
    call PLINK.make_bed_plink1 as extract_chrs{
        input:
            bed_in = normalize_sex_chr_wf.bed_out,
            bim_in = normalize_sex_chr_wf.bim_out,
            fam_in = normalize_sex_chr_wf.fam_out,
            chrs = chrs,
            output_basename = "${output_basename}.extract_chrs",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb
    }

    File norm_bed = extract_chrs.bed_out
    File norm_bim = extract_chrs.bim_out
    File norm_fam = extract_chrs.fam_out

    # Parallelize impute2 id conversion by chr
    scatter(chr_index in range(length(chrs))){
        String chr = chrs[chr_index]

        # Subset bim file to get only bim entries for current chr
        call TSV.tsv_filter as split_bim{
            input:
                tsv_input = norm_bim,
                output_filename = "${output_basename}.chr.${chr}.bim",
                header = false,
                filter_string = "--eq 1:${chr}"
        }

        # Convert IDs to Impute2 format
        call IDCONVERT.convert_variant_ids as impute2_id_bim{
            input:
                in_file = split_bim.tsv_output,
                ref = id_legend_files[chr_index],
                chr = chr,
                in_header = 0,
                in_sep = "tab",
                in_id_col = 1,
                in_chr_col = 0,
                in_pos_col = 3,
                in_a1_col = 4,
                in_a2_col = 5,
                in_missing_allele = in_monomorphic_allele,
                in_deletion_allele = in_deletion_allele,
                ref_deletion_allele = ref_deletion_allele,
                output_filename = "${output_basename}.chr.${chr}.impute2",
                cpu = id_convert_cpu,
                mem_gb = id_convert_mem_gb,
                rescue_rsids = rescue_rsids
        }

        # Mark variants with duplicate IDS if required
        if(remove_duplicates){
            call label_duplicate_variants{
                input:
                    bim_in = impute2_id_bim.output_file,
                    output_basename = "${output_basename}.chr.${chr}.impute2.mrkdup",
                    cpu = duplicate_id_cpu,
                    mem_gb = duplicate_id_mem_gb
            }
        }
        File chr_bim = select_first([label_duplicate_variants.bim_out, impute2_id_bim.output_file])
    }

    # Concatenate chr bims into single bim
    call UTILS.cat as cat_chr_bim{
        input:
            input_files = chr_bim,
            output_filename = "${output_basename}.impute2.bim"
    }

    # Remove duplicates if desired
    if(remove_duplicates){

        # Get list of variant IDs labeled as duplicates
        call get_duplicate_variant_ids{
            input:
                bim_in = cat_chr_bim.output_file,
                output_filename = "${output_basename}.dup_ids.txt"
        }

        # If duplicates, select variant with highest call rate from each group
        if (size(get_duplicate_variant_ids.duplicate_ids) > 0){
            # Get list of variants to remove
            call get_variants_to_remove{
                input:
                    bed_in = norm_bed,
                    bim_in = cat_chr_bim.output_file,
                    fam_in = norm_fam,
                    duplicate_snp_ids = get_duplicate_variant_ids.duplicate_ids,
                    output_basename = output_basename,
                    cpu = plink_cpu,
                    mem_gb = plink_mem_gb
            }

            # Remove duplicates from dataset
            call PLINK.make_bed as remove_dups{
                input:
                    bed_in = norm_bed,
                    bim_in = cat_chr_bim.output_file,
                    fam_in = norm_fam,
                    exclude = get_variants_to_remove.dups_to_remove,
                    output_basename = "${output_basename}.impute2_id.unique",
                    cpu = plink_cpu,
                    mem_gb = plink_mem_gb
            }

            # Remove numbers from duplicate variant IDs
            call fix_ids{
                input:
                    bim_in = remove_dups.bim_out,
                    output_basename = "${output_basename}.impute2_id.unique.fixed"
            }
        }
    }

    output {
        File bed_out = select_first([remove_dups.bed_out, norm_bed])
        File bim_out = select_first([fix_ids.bim_out, cat_chr_bim.output_file])
        File fam_out = select_first([remove_dups.fam_out, norm_fam])
        File? removed_duplicate_ids = get_variants_to_remove.dups_to_remove
    }
}
