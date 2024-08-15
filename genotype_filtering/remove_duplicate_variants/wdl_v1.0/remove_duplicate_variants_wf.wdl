version 1.0

import "plink.wdl" as PLINK

workflow remove_duplicate_variants_wf{

    input {

        File bed_in
        File bim_in
        File fam_in
        String output_basename

        Int label_duplicate_variants_cpu = 1
        Int label_duplicate_variants_mem_gb = 2

        Int plink_cpu = 1
        Int plink_mem_gb = 2

        String image_source = "docker"
        String? ecr_repo

    }

    # Append suffixes to duplicate variants
    call label_duplicate_variants{
        input:
            bim_in = bim_in,
            output_basename = output_basename,
            cpu = label_duplicate_variants_cpu,
            mem_gb = label_duplicate_variants_mem_gb,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    # Get list of variant IDs labeled as duplicates
    call get_duplicate_variant_ids{
        input:
            bim_in = label_duplicate_variants.bim_out,
            output_basename = output_basename,
            image_source = image_source,
            ecr_repo = ecr_repo
    }

    if (size(get_duplicate_variant_ids.duplicate_ids) > 0){
        # Get list of variants to remove
        call get_variants_to_remove{
            input:
                bed_in = bed_in,
                bim_in = label_duplicate_variants.bim_out,
                fam_in = fam_in,
                duplicate_ids = get_duplicate_variant_ids.duplicate_ids,
                output_basename = output_basename,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Remove duplicates from dataset
        call PLINK.make_bed as remove_dups{
            input:
                bed_in = bed_in,
                bim_in = label_duplicate_variants.bim_out,
                fam_in = fam_in,
                exclude = get_variants_to_remove.dups_to_remove,
                output_basename = output_basename,
                cpu = plink_cpu,
                mem_gb = plink_mem_gb,
                image_source = image_source,
                ecr_repo = ecr_repo
        }

        # Remove numbers from duplicate variant IDs
        call fix_ids{
            input:
                bim_in = remove_dups.bim_out,
                output_basename = output_basename,
                image_source = image_source,
                ecr_repo = ecr_repo
        }
    }

    output {
        File bed_out = select_first([remove_dups.bed_out, bed_in])
        File bim_out = select_first([fix_ids.bim_out, bim_in])
        File fam_out = select_first([remove_dups.fam_out, fam_in])
        File? removed_duplicate_ids = get_variants_to_remove.dups_to_remove
    }
}

task label_duplicate_variants{

    input {

        File bim_in
        String output_basename

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 1

    }

    command <<<
        set -e
        # Get sorted list of all variants
        perl -lane 'print $F[1];' ~{bim_in} | sort > tmp.variants.sorted

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
            print join("\t", @F);' ~{bim_in} > ~{output_basename}_duplicates_labeled.bim
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File bim_out = "~{output_basename}_duplicates_labeled.bim"
    }
}

task get_duplicate_variant_ids{

    input {

        File bim_in
        String output_basename

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 1

    }

    command <<<
        grep ___ ~{bim_in} |
            perl -lane 'print $F[1]'> ~{output_basename}_duplicate_ids.txt
        touch ~{output_basename}_duplicate_ids.txt
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File duplicate_ids = "~{output_basename}_duplicate_ids.txt"
    }
}

task get_variants_to_remove{

    input {

        File bed_in
        File bim_in
        File fam_in
        File duplicate_ids
        String? chr
        String output_basename
        String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

        # Runtime environment
        String docker_image = "rtibiocloud/plink:v1.9_2706f72"
        String ecr_image = "rtibiocloud/plink:v1.9_2706f72"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 4
        Int mem_gb = 8

    }

    command <<<
        set -e

        mkdir plink_input

        # Bed file preprocessing
        if [[ ~{bed_in} =~ \.gz$ ]]; then
            unpigz -p ~{cpu} -c ~{bed_in} > plink_input/~{input_prefix}.bed
        else
            # Otherwise just create softlink with normal
            ln -s ~{bed_in} plink_input/~{input_prefix}.bed
        fi

        # Bim file preprocessing
        if [[ ~{bim_in} =~ \.gz$ ]]; then
            unpigz -p ~{cpu} -c ~{bim_in} > plink_input/~{input_prefix}.bim
        else
            ln -s ~{bim_in} plink_input/~{input_prefix}.bim
        fi

        # Fam file preprocessing
        if [[ ~{fam_in} =~ \.gz$ ]]; then
            unpigz -p ~{cpu} -c ~{fam_in} > plink_input/~{input_prefix}.fam
        else
            ln -s ~{fam_in} plink_input/~{input_prefix}.fam
        fi

        # Get call rates of duplicate SNPs
        plink --bfile plink_input/~{input_prefix} \
            --missing \
            --extract ~{duplicate_ids} \
            --out ~{output_basename}_missing

        # Create remove list containing duplicate with higher missing rate
        tail -n +2 ~{output_basename}_missing.lmiss |
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
                }' > ~{output_basename}_duplicates_to_remove.txt
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File dup_missing_report = "~{output_basename}_missing.lmiss"
        File dups_to_remove = "~{output_basename}_duplicates_to_remove.txt"
    }
}

task fix_ids{

    input {

        File bim_in
        String output_basename

        # Runtime environment
        String docker_image = "ubuntu:22.04@sha256:19478ce7fc2ffbce89df29fea5725a8d12e57de52eb9ea570890dc5852aac1ac"
        String ecr_image = "rtibiocloud/ubuntu:22.04_19478ce7fc2ff"
        String? ecr_repo
        String image_source = "docker"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        Int cpu = 1
        Int mem_gb = 1

    }

    command <<<
        set -e
        perl -lne 's/___\d+//; print;' ~{bim_in} > ~{output_basename}.bim
    >>>

    runtime {
        docker: container_image
        cpu: cpu
        memory: "~{mem_gb} GB"
    }

    output {
        File bim_out = "~{output_basename}.bim"
    }
}
