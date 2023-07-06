task locate_high_ld_regions {
    File bimfile

    String docker
    String cpu = 1
    String mem = 2

    command <<<

        awk '
         {if (($1==5 && $4 >= 43964243 && $4 <= 51464243) || 
             ($1==6 && $4 >= 24892021  && $4 <= 33392022) || 
             ($1==8 && $4 >= 7962590   && $4 <= 11962591)  || 
             ($1==11 && $4 >= 45043424 && $4 <= 57243424))
            {print $2}
         }' ${bimfile} > "high_ld_regions.remove"

    >>>

    output {
        File high_ld_regions = "high_ld_regions.remove"
    }

    runtime {
        docker : docker
        cpu : cpu
        mem : "${mem} GB"
    }

    parameter_meta {
        bimfile: "Plink formatted bimfile."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
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


    String docker
    String cpu = 1
    String mem = 2

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
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
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
        mem: "Amount of RAM in GB for the image."
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

    String docker
    String cpu = 1
    String mem = 2

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
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bedfile: "Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Perform LD pruning genotype data that are in high LD regions. These commands produce a pruned subset of markers that are in approximate linkage equilibrium with each other, writing the IDs to plink.prune.in (and the IDs of all excluded variants to plink.prune.out)."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task merge_pruned {

    Array[File] pruned_files

    String docker
    String cpu = 1
    String mem = 2

    command <<<
        cat ${sep=" " pruned_files} > "chr_all_ld_pruned.prune.in"
    >>>

    output {
        File combined_variants = "chr_all_ld_pruned.prune.in"
    }

    runtime {
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
    }

    parameter_meta {
        pruned_files: "Array of chromosome files that contain variants that are in approximate linkage equilibrium (they are independent)."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
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

    String docker
    String cpu = 1
    String mem = 2

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
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bedfile: "Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
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

    String docker
    String cpu = 1
    String mem = 2

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
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        bimfile: "Plink formatted bim file."
        famfile: "Plink formatted fam file."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
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

    String docker
    String cpu
    String mem

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
        File eigenvectors = "${study_name}_${ancestry}_ld_pruned.pca.evec"
    }

    runtime {
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
    }

    parameter_meta {
        bedfile: "LD-pruned, Plink formatted bed file. Not to be confused with UCSC\'s Browser Extensible Data (BED) format."
        bimfile: "LD-pruned, Plink formatted bim file - renamed variant IDs."
        famfile: "LD-pruned, Plink formatted fam file - renamed sample IDs."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Run smartpca on the LD-pruned genotype data."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}

task create_final_file {

    String study_name
    String ancestry

    File evec_file
    File famfile
    String final_file_name = "${study_name}_${ancestry}_ld_pruned_top10_pcs.txt"

    String docker
    String cpu
    String mem

    command <<<

        ## Extract eigenvectors (top 10 PCs) ##
        ## The Principal components are the eigenvectors of the covariance matrix.
        echo "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > top10_pcs.tmp

        tail -n +2 ${evec_file} | \
          awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' >> top10_pcs.tmp

        ## Add original IDs back in
        echo "fid iid" > fid_iid.txt
        awk '{print $1, $2}' ${famfile} >> fid_iid.txt

        paste -d " " fid_iid.txt top10_pcs.tmp > ${final_file_name}
    >>>

    output {
        File final_file = "${final_file_name}"
    }

    runtime {
        docker: docker
        cpu: cpu
        mem: "${mem} GB"
    }

    parameter_meta {
        study_name: "Name of the study/cohort whose genotype data you are analyzing."
        ancestry: "Ancestry of the genotype data you analyzing."
        evec_file: "Top eigenvectors from the PCA."
        famfile: "Original famfile that has the original sample IDs so we can add them back in."
        final_file_name: "Name of the final PCA file that contains the top10 PCs."
        docker: "Docker image"
        cpu: "Number of CPUs for the image."
        mem: "Amount of RAM in GB for the image."
    }

    meta {
        description: "Create top-10 PCs file."
        author: "Jesse Marks"
        email: "jmarks@rti.org"
    }
}
