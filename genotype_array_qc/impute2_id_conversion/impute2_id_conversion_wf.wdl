import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/biocloud_wdl_tools/convert_to_impute2_ids/convert_to_impute2_ids.wdl" as IDCONVERT

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
    String output_basename
    String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

    # Runtime environment
    String docker = "rtibiocloud/plink:v1.9-9e70778"
    Int cpu = 4
    Int mem_gb = 8

    command <<<
        set -e

        mkdir plink_input

        # Bed file preprocessing
        if [[ ${bed_in} =~ \.gz$ ]]; then
            # Append gz tag to let plink know its gzipped input
            gunzip -c ${bed_in} > plink_input/${input_prefix}.bed
        else
            # Otherwise just create softlink with normal
            ln -s ${bed_in} plink_input/${input_prefix}.bed
        fi

        # Bim file preprocessing
        if [[ ${bim_in} =~ \.gz$ ]]; then
            gunzip -c ${bim_in} > plink_input/${input_prefix}.bim
        else
            ln -s ${bim_in} plink_input/${input_prefix}.bim
        fi

        # Fam file preprocessing
        if [[ ${fam_in} =~ \.gz$ ]]; then
            gunzip -c ${fam_in} > plink_input/${input_prefix}.fam
        else
            ln -s ${fam_in} plink_input/${input_prefix}.fam
        fi

        # Get call rates of duplicate SNPs
        plink --bfile plink_input/${input_prefix} \
            --missing \
            --extract ${duplicate_snp_ids} \
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
    String file_in_monomorphic_allele = "0"

    # Resources subtasks
    # The only one you'll likely need to play for huge files is the id_convert cpu/mem
    # But mainly you just need to make sure the mem is roughly the size of the largest id_legend file
    Int split_bed_cpu = 1
    Int split_bed_mem_gb = 2

    Int merge_bed_cpu = 4
    Int merge_bed_mem_gb = 8

    Int id_convert_cpu = 2
    Int id_convert_mem_gb = 6

    Int duplicate_id_cpu = 2
    Int duplicate_id_mem_gb = 6

    Int x_merge_cpu = 1
    Int x_merge_mem_gb = 2

    # Merge X chr to ensure PAR/NONPAR are not split (split-x will fail for pre-split files)
    call PLINK.make_bed as merge_x_chr{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            output_basename = "${output_basename}.mergex",
            merge_x = true,
            merge_no_fail = no_fail,
            cpu = x_merge_cpu,
            mem_gb = x_merge_mem_gb
    }

    # Split X chr by PAR/NONPAR
    call PLINK.make_bed as split_x_chr{
        input:
            bed_in = merge_x_chr.bed_out,
            bim_in = merge_x_chr.bim_out,
            fam_in = merge_x_chr.fam_out,
            output_basename = "${output_basename}.splitx",
            split_x = true,
            build_code = build_code,
            split_no_fail = no_fail,
            cpu = x_merge_cpu,
            mem_gb = x_merge_mem_gb
    }

    # Parallelize impute2 id conversion by chr
    scatter(chr_index in range(length(chrs))){
        String chr = chrs[chr_index]

        # Split by chr
        call PLINK.make_bed as split_bed{
            input:
                bed_in = split_x_chr.bed_out,
                bim_in = split_x_chr.bim_out,
                fam_in = split_x_chr.fam_out,
                output_basename = "${output_basename}.chr.${chr}",
                chr = chr,
                cpu = split_bed_cpu,
                mem_gb = split_bed_mem_gb
        }

        # Convert IDs to Impute2 format
        call IDCONVERT.convert_to_impute2_ids as impute2_id_bim{
            input:
                in_file = split_bed.bim_out,
                legend_file = id_legend_files[chr_index],
                file_in_header = 0,
                id_col = 1,
                chr_col = 0,
                pos_col = 3,
                a1_col = 4,
                a2_col = 5,
                file_in_monomorphic_allele = file_in_monomorphic_allele,
                output_filename = "${output_basename}.impute2",
                cpu = id_convert_cpu,
                mem_gb = id_convert_mem_gb
        }

        # Mark variants with duplicate IDS if required
        if(remove_duplicates){
            call label_duplicate_variants{
                input:
                    bim_in = impute2_id_bim.output_file,
                    output_basename = "${output_basename}.impute2.mrkdup",
                    cpu = duplicate_id_cpu,
                    mem_gb = duplicate_id_mem_gb
            }
        }
        File bim_files = select_first([label_duplicate_variants.bim_out, impute2_id_bim.output_file])
    }

    # Merge into single file
    #Array[File] bim_files = if(remove_duplicates) then select_all(label_duplicate_variants.bim_out) else impute2_id_bim.output_file
    call PLINK.merge_beds{
        input:
            bed_in = split_bed.bed_out,
            bim_in = bim_files,
            fam_in = split_bed.fam_out,
            output_basename = "${output_basename}.impute2_id",
            cpu = merge_bed_cpu,
            mem_gb = merge_bed_mem_gb
    }

    # Remove duplicates if specified
    if(remove_duplicates){
        # Get list of variant IDs labeled as duplicates
        call get_duplicate_variant_ids{
            input:
                bim_in = merge_beds.bim_out,
                output_filename = "${output_basename}.dup_ids.txt"
        }

        # Only do rest of workflow if there are actually duplicate variants to remove
        if (size(get_duplicate_variant_ids.duplicate_ids) > 0){
            # Get list of variants to remove
            call get_variants_to_remove{
                input:
                    bed_in = merge_beds.bed_out,
                    bim_in = merge_beds.bim_out,
                    fam_in = merge_beds.fam_out,
                    duplicate_snp_ids = get_duplicate_variant_ids.duplicate_ids,
                    output_basename = output_basename
            }

            # Remove duplicate SNPs with lowest call rates
            call PLINK.make_bed as remove_dups{
                input:
                    bed_in = merge_beds.bed_out,
                    bim_in = merge_beds.bim_out,
                    fam_in = merge_beds.fam_out,
                    exclude = get_variants_to_remove.dups_to_remove,
                    output_basename = "${output_basename}.impute2_id.unique"
            }

            # Remove numbers from duplicate variant IDs
            call fix_ids{
                input:
                    bim_in = remove_dups.bim_out,
                    output_basename = "${output_basename}.impute2_id.unique"
            }
        }
    }

    output{
        File bed_out = select_first([remove_dups.bed_out, merge_beds.bed_out])
        File bim_out = select_first([fix_ids.bim_out, merge_beds.bim_out])
        File fam_out = select_first([remove_dups.fam_out, merge_beds.fam_out])
        File? all_duplicate_ids = get_duplicate_variant_ids.duplicate_ids
        File? duplicate_id_call_rates = get_variants_to_remove.dup_missing_report
        File? removed_duplicate_ids = get_variants_to_remove.dups_to_remove
    }
}