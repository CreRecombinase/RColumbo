library(gdsfmt)
library(SNPRelate)
library(RColumbo)
library(dplyr)
library(readr)
library(tidyr)

#gdsfile <- "/media/nwknoblauch/Data/GTEx/GTEx_SNP_GDS/Adipose_Subcutaneous.GDS"
#h5file <- "/media/nwknoblauch/Data/GTEx/GTEx_EXP_h5/Adipose_Subcutaneous.h5"
#beddir <- "/media/nwknoblauch/Data/GTEx/GTEx_BED/"
#gctadir <- "/home/nwknoblauch/Downloads/gcta/gcta64"

args <- commandArgs(trailingOnly = T)
gdsfile <- args[1]
bed_pref <-args[2]




#dgn_bed <- "/media/nwknoblauch/Data/DGN/Lev/DGN"
#gtex_h5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data.h5"
#gtex_bed <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed"
#dgn_famf <- "/media/nwknoblauch/Data/DGN/Lev/DGN.fam"
#expf <- "/media/nwknoblauch/Data/DGN/Lev/data_used_for_eqtl_study/cis_data.txt"
#gtex_famf <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed.fam"
#expdir <- "/home/nwknoblauch/Desktop/eQTL/Heritability/expdir/"
cat("Creating file ",bed_pref,"\n")
snpgdsGDS2BED(gdsfile,bed_pref)





