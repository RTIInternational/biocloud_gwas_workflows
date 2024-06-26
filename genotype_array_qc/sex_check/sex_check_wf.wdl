import "biocloud_gwas_workflows/biocloud_wdl_tools/plink/plink.wdl" as PLINK
import "biocloud_gwas_workflows/genotype_array_qc/ld_pruning/ld_prune_wf.wdl" as LD
import "biocloud_gwas_workflows/genotype_array_qc/normalize_sex_chr/normalize_sex_chr_wf.wdl" as NORM
import "biocloud_gwas_workflows/biocloud_wdl_tools/utils/utils.wdl" as UTILS

workflow sex_check_wf{
    File bed_in
    File bim_in
    File fam_in
    File phenotype_file
    String output_basename

    # Phenotype file attributes to be used for re-factoring
    Int header_rows
    Int fid_col
    Int iid_col
    Int sex_col
    String delimiter

    # Args for splitting chrX into PAR/NONPAR
    String build_code
    Boolean no_fail

    # Args for LD pruning
    File? ld_exclude_regions
    String ld_type = "indep-pairphase"
    Int window_size
    Int step_size
    Float r2_threshold
    Float? min_ld_maf
    String? window_size_unit
    Int? x_chr_mode

    # Args for sex check
    Float female_max_f = 0.2
    Float male_min_f = 0.8

    # Runtime options
    Int plink_cpu = 1
    Int plink_mem_gb = 2
    Int ld_cpu = 8
    Int ld_mem_gb = 16
    Int sex_check_cpu = 8
    Int sex_check_mem_gb = 16

    # Container
    String container_source = "docker"

    # Reformat phenotype file
    call format_phenotype_file{
        input:
            phenotype_in = phenotype_file,
            output_filename = "${output_basename}.sexfile.txt",
            header_rows = header_rows,
            fid_col = fid_col,
            iid_col = iid_col,
            sex_col = sex_col,
            delimiter = delimiter,
            container_source = container_source
    }

    # Update sex based on pheno file
    call PLINK.make_bed as update_sex{
        input:
            bed_in = bed_in,
            bim_in = bim_in,
            fam_in = fam_in,
            update_sex = format_phenotype_file.phenotype_out,
            output_basename = "${output_basename}.updated_sex",
            cpu = plink_cpu,
            mem_gb = plink_mem_gb,
            container_source = container_source
    }

    # Make sure PAR/NONPAR regions are split
    call NORM.normalize_sex_chr_wf{
        input:
            bed_in = update_sex.bed_out,
            bim_in = update_sex.bim_out,
            fam_in = update_sex.fam_out,
            expected_chrs = ["23","25"],
            build_code = build_code,
            no_fail = no_fail,
            output_basename = "${output_basename}.x_split",
            plink_cpu = plink_cpu,
            plink_mem_gb = plink_mem_gb,
            container_source = container_source
    }

    # Count number of chrX variants
    call count_chr_x_variants{
        input:
            bim = normalize_sex_chr_wf.bim_out
    }
    Int chr_x_variant_count = count_chr_x_variants.chr_x_variant_count

    if(chr_x_variant_count > 0){
        # Get LD pruned set
        call LD.ld_prune_wf as ld_prune{
            input:
                bed_in = normalize_sex_chr_wf.bed_out,
                bim_in = normalize_sex_chr_wf.bim_out,
                fam_in = normalize_sex_chr_wf.fam_out,
                chr = "23, 25",
                output_basename = "${output_basename}.ldprune",
                ld_type = ld_type,
                window_size = window_size,
                step_size = step_size,
                window_size_unit = window_size_unit,
                r2_threshold = r2_threshold,
                x_chr_mode = x_chr_mode,
                cpu = ld_cpu,
                mem_gb = ld_mem_gb,
                maf = min_ld_maf,
                exclude_regions = ld_exclude_regions,
                container_source = container_source
        }

        # Do Sex Check
        call PLINK.sex_check{
            input:
                bed_in = ld_prune.bed_out,
                bim_in = ld_prune.bim_out,
                fam_in = ld_prune.fam_out,
                female_max_f = female_max_f,
                male_min_f = male_min_f,
                output_basename = "${output_basename}.sexcheck",
                cpu = sex_check_cpu,
                mem_gb = sex_check_mem_gb,
                container_source = container_source
        }
    }

    if(chr_x_variant_count == 0){
        call create_empty_file as create_output_file {
            input:
                file_name = "${output_basename}.no_sexcheck.output",
                container_source = container_source
        }

        call create_empty_file as create_problems_file {
            input:
                file_name = "${output_basename}.no_sexcheck.problems",
                container_source = container_source
        }

        call create_empty_file as create_remove_file {
            input:
                file_name = "${output_basename}.no_sexcheck.remove",
                container_source = container_source
        }
    }

    File sex_check_output = select_first([sex_check.plink_sex_check_output, create_output_file.output_file])
    File problems = select_first([sex_check.sex_check_problems, create_problems_file.output_file])
    File remove = select_first([sex_check.samples_to_remove, create_remove_file.output_file])

    # Output sex check output and list of samples with sex mismatches between phenotype file and sex check
    output{
        File plink_sex_check_output = sex_check_output
        File sex_check_problems = problems
        File samples_to_remove = remove
    }
}
task format_phenotype_file{
    File phenotype_in
    String output_filename
    Int header_rows
    Int fid_col
    Int iid_col
    Int sex_col
    String delimiter

    Int tail_n = header_rows + 1

    # Runtime environment
    String docker = "ubuntu:22.04@sha256:a6d2b38300ce017add71440577d5b0a90460d0e57fd7aec21dd0d1b0761bbfb2"
    String ecr = "public.ecr.aws/ubuntu/ubuntu:22.04_stable"
    String container_source = "docker"
    String container_image = if(container_source == "docker") then docker else ecr
    Int cpu = 1
    Int mem_gb = 1

    command {
        set -e

        if [[ ${phenotype_in} =~ \.gz$ ]]; then
            gunzip -c ${phenotype_in} > ${output_filename}.tmp
        else
            ln -s ${phenotype_in} ${output_filename}.tmp
        fi

        tail -n +${tail_n} ${output_filename}.tmp |
        perl -lne '
            $delimiter = lc("${delimiter}");
            $delimiter = ($delimiter eq "comma") ? "," : (($delimiter eq "tab") ? "\t" : (($delimiter eq "space") ? " " : ""));
            chomp;
            @F = split($delimiter);
            print join("\t", $F[${fid_col}], $F[${iid_col}], $F[${sex_col}])' > ${output_filename}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File phenotype_out = "${output_filename}"
    }
}

task count_chr_x_variants{
    File bim

    # Runtime environment
    String docker = "ubuntu:22.04@sha256:a6d2b38300ce017add71440577d5b0a90460d0e57fd7aec21dd0d1b0761bbfb2"
    String ecr = "public.ecr.aws/ubuntu/ubuntu:22.04_stable"
    String container_source = "docker"
    String container_image = if(container_source == "docker") then docker else ecr
    Int cpu = 1
    Int mem_gb = 1

    command{
        grep -P "^(23|25)" ${bim} | wc -l | cut -d" " -f1
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Int chr_x_variant_count = read_int(stdout())
    }
}

task create_empty_file{
    String file_name

    # Runtime environment
    String docker = "ubuntu:22.04@sha256:a6d2b38300ce017add71440577d5b0a90460d0e57fd7aec21dd0d1b0761bbfb2"
    String ecr = "public.ecr.aws/ubuntu/ubuntu:22.04_stable"
    String container_source = "docker"
    String container_image = if(container_source == "docker") then docker else ecr
    Int cpu = 1
    Int mem_gb = 1

    command {
        set -e
        touch ${file_name}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File output_file = "${file_name}"
    }
}

