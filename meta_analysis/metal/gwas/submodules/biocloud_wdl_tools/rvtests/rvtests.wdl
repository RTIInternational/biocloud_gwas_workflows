task rvtests {
    File inVCF
    File phenoFile
    String output_basename
    File? covarFile
    Boolean? inverseNormal
    Boolean? useResidualAsPhenotype
    Boolean? outputRaw
    Boolean? sex
    Boolean? qtl
    Boolean? multipleAllele
    Boolean? xHemi
    String? xLabel
    String? xParRegion
    String? dosage
    String? phenoName
    Int? mpheno

    Array[String]? covarsMaybe
    Array[String]? singleTestsMaybe
    Array[String]? burdenTestsMaybe
    Array[String]? vtTestsMaybe
    Array[String]? kernelTestsMaybe
    Array[String]? metaTestsMaybe

    String covarsPrefix = if defined(covarsMaybe) then "--covar-name "  else ""
    String singlePrefix = if defined(singleTestsMaybe) then "--single " else ""
    String burdenPrefix = if defined(burdenTestsMaybe) then "--burden " else ""
    String vtPrefix = if defined(vtTestsMaybe) then "--vt " else ""
    String kernelPrefix =if defined(kernelTestsMaybe) then "--kernel " else ""
    String metaPrefix = if defined(metaTestsMaybe) then "--meta " else ""

    File? peopleIncludeFile
    Float? freqUpper
    Float? freqLower
    File? rangeFile
    File? siteFile
    String? site
    Int? siteMACMin
    Int ?siteDepthMin
    Int ?siteDepthMax
    String? annoType
    String? impute
    Boolean? imputePheno
    File? geneFile
    Array[String]? genes
    String genesPrefix = if defined(genes) then "--genes " else ""
    File? setFile
    Array[String]? set
    String setPrefix = if defined(set) then "--set " else ""

    File? kinship
    File? xHemiKinship
    File? kinshipEigen
    File? xHemiKinshipEigen
    Boolean? hideCovar
    Boolean? outputID

    # Runtime attributes
    String docker = "rtibiocloud/rvtests:v2.1.0-8d966cb"
    Int cpu = 1
    Int mem_gb = 2
    Int max_retries = 3

    command {
        rvtest --inVcf ${inVCF} \
	        --out ${output_basename} \
	        --pheno ${phenoFile} \
		    ${"--pheno-name " + phenoName} \
		    ${"--mpheno " + mpheno} \
		    ${"--covar " + covarFile } \
		    ${"--xLabel " + xLabel } \
		    ${"--xParRegion " + xParRegion } \
		    ${"--dosage " + dosage } \
		    ${true="--outputRaw" false="" outputRaw} \
		    ${true="--sex" false="" sex} \
		    ${true="--qtl" false="" qtl} \
		    ${true="--multipleAllele" false="" multipleAllele} \
		    ${true="--inverseNormal" false="" inverseNormal} \
		    ${true="--useResidualAsPhenotype" false="" useResidualAsPhenotype} \
            ${true="--xHemi" false="" xHemi }\
            ${covarsPrefix} ${sep="," covarsMaybe} \
            ${singlePrefix} ${sep="," singleTestsMaybe} \
            ${burdenPrefix} ${sep="," burdenTestsMaybe} \
            ${vtPrefix} ${sep="," vtTestsMaybe} \
            ${kernelPrefix} ${sep="," kernelTestsMaybe} \
            ${metaPrefix} ${sep="," metaTestsMaybe} \
            ${ "--peopleIncludeFile " + peopleIncludeFile } \
            ${ "--freqUpper " + freqUpper} \
            ${ "--freqLower " + freqLower} \
            ${ "--rangeFile " + rangeFile } \
            ${ "--siteFile " +  siteFile } \
            ${ "--siteMACMin " +  siteMACMin } \
            ${ "--siteDepthMin " +  siteDepthMin } \
            ${ "--siteDepthMax " +  siteDepthMax } \
            ${ "--annoType " +  annoType } \
            ${ "--impute " +  impute } \
            ${true="--imputePheno" false="" imputePheno} \
            ${ "--geneFile " +  geneFile } \
            ${genesPrefix} ${sep="," genes} \
            ${ "--setFile " +  setFile } \
            ${setPrefix} ${sep="," set} \
            ${ "--kinship " +  kinship } \
            ${ "--xHemiKinship " +  xHemiKinship } \
            ${ "--kinshipEigen " +  kinshipEigen } \
            ${ "--xHemiKinshipEigen " +  xHemiKinshipEigen } \
            ${true="--hide-covar" false="" hideCovar} \
		    ${true="--outputID" false="" outputID} \
            ${"--numThread " + cpu }
    }
    output {
        Array[File] assoc_outputs = glob( "${output_basename}*.assoc*.gz" )
        File log_file = "${output_basename}.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

}

task vcf2kinship  {
    File? inputVcf
    File? pedfile
    String? dosage
    Boolean? xHemi
    String? xLabel
    Float? maxMiss
    Float? minMAF
    Float? minSiteQual

    String output_basename
    Boolean useBaldingNicols
    Boolean useIBS

    # Runtime attributes
    String docker = "rtibiocloud/rvtests:v2.1.0-8d966cb"
    Int cpu = 4
    Int mem_gb = 8
    Int max_retries = 3

    command  {
        vcf2kinship  ${"--inVcf " + inputVcf} \
            ${"--ped " + pedfile} \
            ${true="--bn" false="" useBaldingNicols} \
            ${true="--ibs" false=""  useIBS} \
            ${true="--xHemi" false="" xHemi} \
            ${"--dosage " + dosage } \
            ${"--xLabel " + xLabel } \
            ${"--maxMiss " + maxMiss } \
            ${"--minMAF " + minMAF } \
            ${"--minSiteQual " + minSiteQual } \
            --thread ${cpu} \
            --out ${output_basename}

        # Hack because WDL doesn't allow optional output files
        touch ${output_basename}.kinship
        touch ${output_basename}.xHemi.kinship
    }

    output {
        File kinship_matrix = "${output_basename}.kinship"
        File xHemi_kinship_matrix = "${output_basename}.xHemi.kinship"
        File kinship_log = "${output_basename}.vcf2kinship.log"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }
}

task combineKinship  {
    Array[File] kinship_matrices
    Array[File] vcf2kinship_logs
    String output_basename

    # Runtime attributes
    String docker = "rtibiocloud/rvtests:v2.1.0-8d966cb"
    Int cpu = 16
    Int mem_gb = 16
    Int max_retries = 3

    command  {
        # Apparently the logs have to be in the same damn directory as the kinship mats

        # Copy kinship mats to working directory
        for file in ${sep=" " kinship_matrices}; do
            cp $file .
        done

        # Copy log files to working directory
        for file in ${sep=" " vcf2kinship_logs}; do
            cp $file .
        done

        combineKinship \
            --out ${output_basename} \
            --thread ${cpu} \
            ./*.kinship
    }

    output {
        File kinship_matrix = "${output_basename}.kinship"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }
}
