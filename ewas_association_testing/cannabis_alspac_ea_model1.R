#!/usr/bin/Rscript

args <- commandArgs(TRUE)

loop = TRUE
target_filenames <- list()

# read in the phenotype file, DNAm data, and output directory.
while (loop) {
    
        if (args[1] == "--pheno") {
            pheno_file = args[2]
        }

        if (args[1] == "--sample_name") {
            sample_name = args[2]
        }

        if (args[1] == "--methRdata") {
            methRdata = args[2]
        }

        if (args[1] == "--output") {
            output = args[2]
        }

        if (length(args)>1) {
            args<-args[2:length(args)]
        } else{
              loop = FALSE
        }
} 


# these libraries have already been installed on the EWAS docker image.
library(parallel)  # to use multicore approach - part of base R – can be omitted if lapply() used instead of mclapply()
library(MASS) # rlm function for robust linear regression
library(lmtest) #to use coeftest
library(sandwich) #Huberís estimation of the standard error
library(data.table)
library(stats)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 

####################################################################################################

LMtest = function(meth_matrix, methcol, outcome, X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13) {
        mod = try(lm(meth_matrix[, methcol]~outcome+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13))
	#right_side <- paste("outcome", "X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13", collapse="+")
	#print(right_side)
        #mod = try(lm(meth_matrix[, methcol]~right_side))

	if(class(mod) == "try-error"){
	print(paste("error thrown by column", methcol))
	invisible(rep(NA, 3)) # doesn't print if not assigned. Use in place of return

	}else cf = coef(summary(mod))
	cf[2, c("Estimate", "Std. Error", "Pr(>|t|)")] 
	}

# read in phenotype file
cat("Reading phenotype data......\n")
pheno <- read.csv(pheno_file, header=T, stringsAsFactors=F, sep=" ")
# remove NA from phenotypes
pheno[is.na(pheno$cannabisUse),]$cannabisUse <- 0

pheno4EWAS<- pheno
cat("Phenotype data has ",dim(pheno)[1]," rows and ",dim(pheno)[2]," columns.\n\n")

cat("Loading DNA methylation data......\n")
    load(methRdata)
    
beta_matrix <- t(bVals_chr[,colnames(bVals_chr) %in% pheno4EWAS[, sample_name]])

cat("DNAm data has ",dim(beta_matrix)[1]," rows and ",dim(beta_matrix)[2]," columns.\n\n")

    #Run adjusted EWAS
    ##modify with covariate names
    #  Sample_Name cidB3176 ALN qlet age_at_DNAm cannabisUse smoking drinking BMI Bcell CD4T CD8T Gran Mono NK sv1 sv2 sv3 sv4 sv5 sv6
    system.time(ind.res <- mclapply(setNames(seq_len(ncol(beta_matrix)), 
                                         dimnames(beta_matrix)[[2]]), 
                                LMtest, 
                                meth_matrix=beta_matrix, 
                                outcome=pheno4EWAS$cannabisUse, X1=pheno4EWAS$age_at_DNAm, 
                                X2=pheno4EWAS$sv1, X3=pheno4EWAS$sv2, 
                                X4=pheno4EWAS$sv3, X5=pheno4EWAS$sv4, 
                                X6=pheno4EWAS$sv5, X7=pheno4EWAS$sv6,
                                X8=pheno4EWAS$Bcell, X9=pheno4EWAS$CD4T, 
                                X10=pheno4EWAS$CD8T, X11=pheno4EWAS$Mono, 
                                X12=pheno4EWAS$Gran, X13=pheno4EWAS$NK)) 


setattr(ind.res, 'class', 'data.frame')
setattr(ind.res, "row.names", c(NA_integer_,4))
setattr(ind.res, "names", make.names(names(ind.res), unique=TRUE))
probelistnames <- names(ind.res)
all.results <- t(data.table(ind.res))
all.results<-data.table(all.results)
all.results[, probeID := probelistnames]
setnames(all.results, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(all.results, c("probeID","BETA","SE", "P_VAL"))
rm(probelistnames, ind.res)

#Add column for number of samples for each probe
tbeta_matrix<-t(beta_matrix) #transform methylation data again so that rows are probes and columns are samples
all.results<-all.results[match(rownames(tbeta_matrix),all.results$probeID),] # match order of all.results with order of probes in tbeta_matrix
all.results$N<- rowSums(!is.na(tbeta_matrix))

#Add columns to include Illumina probe annotation
annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
a<-as.data.frame(all.results)
b<-as.data.frame(annEPIC)
all.results.annot<-merge(a, b, by.x="probeID", by.y="Name")
all.results.annot<-all.results.annot[order(all.results.annot$P_VAL),] #sort by P_VAL
write.table(all.results.annot, 
	    paste(output, "_",all.results.annot$chr[1],"_",Sys.Date(),".csv", sep = ""), #paste(out_dir,"Cannabis_ALSPAC_Model1_ANNOT_RESULTS_",all.results.annot$chr[1],"_",Sys.Date(),".csv", sep = ""),
            na="NA",sep = ",",row.names=FALSE) #Export full results; these will be used later for plotting
