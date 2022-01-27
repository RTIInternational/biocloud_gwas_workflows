import "biocloud_wdl_tools/gem/gem.wdl" as GEM

workflow test_gem{

    # GEM Input/Output File Options:
    File pheno_file
    String out
    File? bgen
    File? sample
    File? pgen
    File? pvar
    File? psam
    File? bed
    File? bim
    File? fam
    String? output_style

    # GEM Phenotype File Options:
    String sampleid_name
    String pheno_name
    String? exposure_names
    String? int_covar_names
    String? covar_names
    Int? robust
    Float? tol
    String? delim
    String? missing_value
    Int? center
    Int? scale
    String? categorical_names
    Int? cat_threshold

    # GEM Filtering Options:
    Float? maf
    Float? miss_geno_cutoff
    String? include_snp_file

    # GEM Performance Options:
    Int? threads
    Int? stream_snps

    # GEM Runtime attributes
    Int gem_cpu = 1
    Int gem_mem_gb = 2

    call GEM.gem as gem {
        input:
            pheno_file = pheno_file,
            out = out,
            bgen = bgen,
            sample = sample,
            pgen = pgen,
            pvar = pvar,
            psam = psam,
            bed = bed,
            bim = bim,
            fam = fam,
            output_style = output_style,
            sampleid_name = sampleid_name,
            pheno_name = pheno_name,
            exposure_names = exposure_names,
            int_covar_names = int_covar_names,
            covar_names = covar_names,
            robust = robust,
            tol = tol,
            delim = delim,
            missing_value = missing_value,
            center = center,
            scale = scale,
            categorical_names = categorical_names,
            cat_threshold = cat_threshold,
            maf = maf,
            miss_geno_cutoff = miss_geno_cutoff,
            include_snp_file = include_snp_file,
            threads = threads,
            stream_snps = stream_snps,
            cpu = gem_cpu,
            mem_gb = gem_mem_gb
    }

    output{
        File summary_stats = gem.summary_stats
        File log = gem.log
    }
}

