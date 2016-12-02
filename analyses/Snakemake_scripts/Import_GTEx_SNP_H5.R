library(RColumbo)
library(dplyr)
library(methods)
library(rhdf5)

args  <- commandArgs(trailingOnly = T)
#snpgzfile <- "/home/nwknoblauch/Desktop/eQTL/Snake/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices/Adipose_Subcutaneous_Analysis.snps.txt.gz"
#expgzfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_expression_matrices/Adipose_Subcutaneous_Analysis.v6p.normalized.expression.bed.gz"
#covgzfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_covariates/Adipose_Subcutaneous_Analysis.v6p.covariates.txt"
#h5file <- "/media/nwknoblauch/Data/GTEx/GTEx_SNP_h5/Adipose_Subcutaneous.h5"
GDSfile <- args[1]
h5file <- args[2]

gtissue <- gsub(".+/([^/]+).GDS","\\1",GDSfile)
htissue <- gsub(".+/([^/]+).h5","\\1",h5file)
stopifnot(gtissue==htissue)

gdsn_to_h5(geno_gdsn = GDSfile,geno_h5 = h5file,chunksize = 650000,doFlip = F)

