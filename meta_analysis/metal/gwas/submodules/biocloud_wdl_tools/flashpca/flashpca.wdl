task flashpca{
    File bed_in
    File bim_in
    File fam_in
    String input_prefix = basename(sub(bed_in, "\\.gz$", ""), ".bed")

    Int? ndim
    String? standx    # Standardization for genotypes (binom2 | binom)
    String? standy   # Standardization for phenotypes (sd | binom2)
    Int? div         # Whether to divide eigenvalues by p, n-1, or none (p|n1|none)
    Int? tol         # Tolerance for PCA iterations

    Boolean? batch
    Int? blocksize
    Int? seed
    File? pheno
    Int? precision

    Boolean? project
    File? inload
    File? inmaf
    File? inmeansd

    Boolean? verbose
    Boolean? notime
    Boolean? check

    # Runtime environment
    String docker = "rtibiocloud/flashpca:v2.0-9b4c1b9"
    Int cpu = 1
    Int mem_gb = 1
    Int mem_mb = mem_gb * 1000

    command {
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

        # Run flashPCA
        flashpca --bfile plink_input/${input_prefix} \
            --memory ${mem_mb} \
            --numthreads ${cpu} \
            ${'--blocksize ' + blocksize} \
            ${'--seed ' + seed} \
            ${true='--batch' false="" batch} \
            ${'--pheno ' + pheno} \
            ${'--ndim ' + ndim} \
            ${'--standx ' + standx} \
            ${'--standy ' + standy} \
            ${'--div ' + div} \
            ${'--tol ' + tol} \
            ${true='--project' false="" project} \
            ${'--inload ' + inload} \
            ${'--inmaf ' + inmaf} \
            ${'--inmeansd ' + inmeansd} \
            ${'--inload ' + inload} \
            ${'--precision ' + precision} \
            ${true='--verbose' false="" verbose} \
            ${true='--notime' false="" notime} \
            ${true='--check' false="" check} \
            --outload loadings.txt \
            --outmeansd meansd.txt
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File eigenvectors = "eigenvectors.txt"
        File pcs = "pcs.txt"
        File eigenvalues = "eigenvalues.txt"
        File pve = "pve.txt"
        File loadings = "loadings.txt"
        File meansd = "meansd.txt"
    }
}

task get_non_ancestry_informative_snps{
    File pca_loadings
    String output_filename
    Float loading_value_cutoff = 0.003
    Float cutoff_step_size = 0.001
    Float max_cutoff = 0.01
    Int max_snps = 100000
    Int min_snps = 10000

    # Runtime environment
    String docker = "rtibiocloud/tsv-utils:v1.4.4-8d966cb"
    Int cpu = 1
    Int mem_gb = 2

    command<<<
        set -e

        # Initialize empty SNPs file
        touch snps.txt

        # Loading value cutoff for defining 'ancestral informative SNP'
        max_load=${loading_value_cutoff}

        # Incrementally loosen criteria until enough SNPs are found or you go above max_cutoff
        # Awk command used to check max_cutoff is because bash doesn't natively do floating point arithmatic/comparisons
        while [ $(wc -l snps.txt | cut -d" " -f1) -lt ${min_snps} ] && [ $(awk 'BEGIN {print ("'$max_load'" <= ${max_cutoff})}') -eq 1 ]
        do

            echo "Using loading value cutoff of $max_load until unless fewer than ${min_snps} SNPs found!"

            # Filter out any snps with PC loading value > cutoff for any of first 3 PCs
            tsv-filter \
		        --header \
		        --le "3:$max_load" \
		        --le "4:$max_load" \
		        --le "5:$max_load" \
                ${pca_loadings} > snps.txt

            echo "Found $(wc -l snps.txt | cut -d' ' -f1) SNPs with cutoff $max_load"

            # Loosen threshold by step size (using awk bc we don't have 'bc' on this docker image; ugly but oh well)
            max_load=$(awk -v var=$max_load 'BEGIN{print var + ${cutoff_step_size}}')
        done

        # Subset SNPs to least informative if greater than max
        if [ $(wc -l snps.txt | cut -d" " -f1) -gt ${max_snps} ]
        then
            # Sort by averge loading and take top max_snps SNPs
            awk '{if(NR > 1){sum = 0; for (i = 3; i <= NF; i++) sum += $i; sum /= (NF-2); print $1,sum}}' snps.txt | \
            sort -k2 | \
            head -${max_snps} > tmp.txt
            mv tmp.txt snps.txt
        fi

        # Output only SNP ids of ancestrally uninformative SNPs
        tail -n +2 snps.txt | cut -f1 > ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File snps_to_keep = "${output_filename}"
    }

}