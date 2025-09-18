version 1.1

workflow eigenstrat_smartpca{

    input {

        String study_name
        String ancestry

        File bedfile
        File bimfile
        File famfile

        File reference_file

        String image_source = "docker"
        String? ecr_repo

        String cpu = 1
        String mem_gb = 2

    }

    call locate_high_ld_regions {
        input:
            bimfile = bimfile,
            reference_file = reference_file,
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call remove_high_ld_regions {
        input:
            study_name = study_name,
            ancestry = ancestry,
            bedfile = bedfile, 
            bimfile = bimfile,
            famfile = famfile, 
            high_ld_regions = locate_high_ld_regions.snps_in_high_ld_regions,
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call ld_pruning {
        input:
            study_name = study_name,
            ancestry = ancestry,
            bedfile = remove_high_ld_regions.output_bed, 
            bimfile = remove_high_ld_regions.output_bim,
            famfile = remove_high_ld_regions.output_fam, 
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call merge_pruned {
        input:
            pruned_files = ld_pruning.pruned_files,
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call extract_ld_variants {
        input:
            study_name = study_name,
            ancestry = ancestry,
            bedfile = remove_high_ld_regions.output_bed,
            bimfile = remove_high_ld_regions.output_bim,
            famfile = remove_high_ld_regions.output_fam,
            combined_variants = merge_pruned.combined_variants,
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call rename_bimfam {
        input:
            bimfile = extract_ld_variants.output_bim,
            famfile = extract_ld_variants.output_fam,
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    call run_smartpca {
        input:
            ancestry = ancestry,
            study_name = study_name,
            bedfile = extract_ld_variants.output_bed,
            bimfile = rename_bimfam.ids_renamed_bim,
            famfile = rename_bimfam.ids_renamed_fam,
            image_source = image_source,
            ecr_repo = ecr_repo,
            cpu = cpu,
            mem_gb = mem_gb
    }

    output {
        File evec = run_smartpca.evec
        File eval = run_smartpca.eval
        File snpweight = run_smartpca.snpweight
        File log = run_smartpca.log
    }
}

task locate_high_ld_regions {

    input {

        File bimfile
        File reference_file

        String image_source = "docker"
        String? ecr_repo
        String docker_image = "python:3.10"
        String ecr_image = "python:3.10"
        String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"
        
        String cpu = 1
        String mem_gb = 2

    }

    command <<<

      python - <<EOF
      reference_file = "${reference_file}"
      bim_file = "${bimfile}"
      output_file = "snps_to_remove_in_high_ld_regions.txt"

      regions = []
      with open(reference_file) as ref_file:
          next(ref_file)  # Skip header
          for line in ref_file:
              chrom, start, end = line.strip().split()
              regions.append((chrom, int(start), int(end)))

      with open(bim_file) as bim_file, open(output_file, 'w') as out_file:
          for line in bim_file:
              columns = line.strip().split()
              chrom = columns[0]
              snp = columns[1]
              pos = int(columns[3])

              # Check each region to find a match
              for region in regions:
                  if chrom == region[0]:
                      if region[1] < pos <= region[2]:
                          # Write SNP to the output file
                          out_file.write(snp + "\n")
      EOF
    >>>

    output {
        File snps_in_high_ld_regions = "snps_to_remove_in_high_ld_regions.txt"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        bimfile: "Plink formatted bimfile."
        reference_file: "Tab Separated text file containing regions of high LD. Should contain 3 columns and a header: chromosome, start position, and stop position. See https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD) for examples."
        docker: "Docker image with python3, such as python:3.10"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Create a list of SNPs to remove that are in high LD regions."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task remove_high_ld_regions {

    String study_name
    String ancestry

    File bedfile
    File bimfile
    File famfile
    File high_ld_regions

    String image_source = "docker"
    String? ecr_repo
    String docker_image = "rtibiocloud/plink:v1.9_178bb91"
    String ecr_image = "rtibiocloud/plink:v1.9_178bb91"
    String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"

    String cpu = 1
    String mem_gb = 2

    command <<<
        plink \
            --noweb \
            --bed ${bedfile} \
            --bim ${bimfile} \
            --fam ${famfile} \
            --exclude ${high_ld_regions} \
            --make-bed \
            --out ${study_name}_${ancestry}_genotypes_high_ld_regions_removed
    >>>

    output {
        File output_bed = "${study_name}_${ancestry}_genotypes_high_ld_regions_removed.bed"
        File output_bim = "${study_name}_${ancestry}_genotypes_high_ld_regions_removed.bim"
        File output_fam = "${study_name}_${ancestry}_genotypes_high_ld_regions_removed.fam"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bedfile: "Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        high_ld_regions: "List of SNPs to remove that are in high LD regions."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Remove SNPs from genotype data that are in high LD regions."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task ld_pruning {

    String study_name
    String ancestry

    File bedfile
    File bimfile
    File famfile

    String image_source = "docker"
    String? ecr_repo
    String docker_image = "rtibiocloud/plink:v1.9_178bb91"
    String ecr_image = "rtibiocloud/plink:v1.9_178bb91"
    String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"

    String cpu = 1
    String mem_gb = 2

    command <<<
        # --indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
        for chr in {1..22}; do
            plink \
                --noweb \
                --bed ${bedfile} \
                --bim ${bimfile} \
                --fam ${famfile} \
                --indep-pairwise 1500 150 0.2 \
                --chr $chr \
                --out ${study_name}_${ancestry}_genotypes_high_ld_regions_removed_ld_pruned_chr$chr
        done
    >>>

    output {
        Array[File] pruned_files = glob("*.prune.in")
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bedfile: "Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Perform LD pruning genotype data that are in high LD regions. These commands produce a pruned subset of markers that are in approximate linkage equilibrium with each other, writing the IDs to plink.prune.in (and the IDs of all excluded variants to plink.prune.out)."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task merge_pruned {

    Array[File] pruned_files

    String image_source = "docker"
    String? ecr_repo
    String docker_image = "ubuntu:22.04"
    String ecr_image = "ubuntu:22.04"
    String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"

    String cpu = 1
    String mem_gb = 2

    command <<<
        cat ${sep=" " pruned_files} > "chr_all_ld_pruned.prune.in"
    >>>

    output {
        File combined_variants = "chr_all_ld_pruned.prune.in"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        pruned_files: "Array of chromosome files that contain variants that are in approximate linkage equilibrium (they are independent)."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Merge the prune.in files - which are the pruned subset of markers that are in approximate linkage equilibrium with each other."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task extract_ld_variants {

    String study_name
    String ancestry

    File bedfile
    File bimfile
    File famfile

    File combined_variants

    String image_source = "docker"
    String? ecr_repo
    String docker_image = "rtibiocloud/plink:v1.9_178bb91"
    String ecr_image = "rtibiocloud/plink:v1.9_178bb91"
    String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"

    String cpu = 1
    String mem_gb = 2

    command <<<
        plink \
            --noweb \
            --bed ${bedfile} \
            --bim ${bimfile} \
            --fam ${famfile} \
            --extract ${combined_variants} \
            --make-bed \
            --out ${study_name}_${ancestry}_genotypes_high_ld_regions_removed_ld_pruned
    >>>

    output {
        File output_bed = "${study_name}_${ancestry}_genotypes_high_ld_regions_removed_ld_pruned.bed"
        File output_bim = "${study_name}_${ancestry}_genotypes_high_ld_regions_removed_ld_pruned.bim"
        File output_fam = "${study_name}_${ancestry}_genotypes_high_ld_regions_removed_ld_pruned.fam"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bedfile: "Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Perform LD pruning genotype data that are in high LD regions. These commands produce a pruned subset of markers that are in approximate linkage equilibrium with each other, writing the IDs to plink.prune.in (and the IDs of all excluded variants to plink.prune.out)."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task rename_bimfam {

    File bimfile
    File famfile

    String image_source = "docker"
    String? ecr_repo
    String docker_image = "rtibiocloud/plink:v1.9_178bb91"
    String ecr_image = "rtibiocloud/plink:v1.9_178bb91"
    String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"

    String cpu = 1
    String mem_gb = 2

    command <<<
        ## Rename BIM/FAM file IDs.
        awk '{$2="ID_"NR; print $0}' ${bimfile} > "ids_renamed.bim"
        awk '{$1="ID_"NR; $2="ID_"NR; print $0}' ${famfile} > "ids_renamed.fam"
    >>>

    output {
        File ids_renamed_bim = "ids_renamed.bim"
        File ids_renamed_fam = "ids_renamed.fam"
    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Rename the variant IDs in the bim file and the sample IDs in the fam file. smartpca will throw an error if the IDs are too long. So we changed them to a generic ID for the smartpca analysis, then incorporate them back in later."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task run_smartpca {
    String ancestry
    String study_name

    File bedfile
    File bimfile
    File famfile

    String image_source = "docker"
    String? ecr_repo
    String docker_image = "rtibiocloud/eigensoft:v6.1.4_2d0f99b"
    String ecr_image = "rtibiocloud/eigensoft:v6.1.4_2d0f99b"
    String container_image = if(image_source == "docker") then docker_image else "~{ecr_repo}/~{ecr_image}"

    String cpu
    String mem_gb

    command <<<

        /opt/EIG-6.1.4/bin/smartpca.perl \
            -i ${bedfile} \
            -a ${bimfile} \
            -b ${famfile} \
            -o ${study_name}_${ancestry}_ld_pruned.pca \
            -p ${study_name}_${ancestry}_ld_pruned.plot \
            -e ${study_name}_${ancestry}_ld_pruned.eval \
            -l ${study_name}_${ancestry}_ld_pruned.pca.log \
            -m 0
    >>>

    output {

        File evec = "${study_name}_${ancestry}_ld_pruned.pca.evec"
        File eval = "${study_name}_${ancestry}_ld_pruned.eval"
        File snpweight = "${study_name}_${ancestry}_ld_pruned.snpweight"
        File log = "${study_name}_${ancestry}_ld_pruned.pca.log"

    }

    runtime {
        docker: container_image
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        bedfile: "LD-pruned, Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "LD-pruned, Plink formatted bim file - renamed variant IDs."
        famfile: "LD-pruned, Plink formatted fam file - renamed sample IDs."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem_gb: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Run smartpca on the LD-pruned genotype data."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}
