library(RColumbo)
library(dplyr)
library(methods)
library(rhdf5)

args  <- commandArgs(trailingOnly = T)
#snpgzfile <- "/home/nwknoblauch/Desktop/eQTL/Snake/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices/Adipose_Subcutaneous_Analysis.snps.txt.gz"
#expgzfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_expression_matrices/Adipose_Subcutaneous_Analysis.v6p.normalized.expression.bed.gz"
#covgzfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_covariates/Adipose_Subcutaneous_Analysis.v6p.covariates.txt"
#h5file <- "/media/nwknoblauch/Data/GTEx/GTEx_SNP_h5/Adipose_Subcutaneous.h5"
snpgzfile <- args[1]
expgzfile <- args[2]
covgzfile <- args[3]
h5file <- args[4]

stissue <- gsub(".+/([^/]+)_Analysis.snps.txt.gz","\\1",snpgzfile)
etissue <- gsub(".+/([^/]+)_Analysis.v6p.normalized.expression.bed.gz","\\1",expgzfile)
ctissue <- gsub(".+/([^/]+)_Analysis.v6p.covariates.txt","\\1",covgzfile)
htissue <- gsub(".+/([^/]+).h5","\\1",h5file)
stopifnot(stissue==etissue,etissue==ctissue,ctissue==htissue)


gc()
#colnum <- length(scan(snpgzfile,what=character(),sep="\t",nlines = 1)[-1])
if(file.exists(h5file)){
  file.remove(h5file)
}
#write_genotype_h5(snpgzfile,Nind = colnum,chunksize = 650000,h5file = h5file,doFlip = F,deflate_level = 4)
gc()
write_covar_h5(covgzfile,h5file,chunksize=1,deflate_level=4)
gc()
read_gtex_expression(expgzfile,h5file)
gc()

