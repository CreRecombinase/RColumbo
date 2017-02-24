library(gdsfmt)
library(SNPRelate)
library(RColumbo)
library(dplyr)
library(readr)
library(tidyr)

# gdsfile <- "/media/nwknoblauch/Data/GTEx/GTEx_SNP_GDS/Adipose_Subcutaneous.GDS"
# h5file <- "/media/nwknoblauch/Data/GTEx/GTEx_EXP_h5/Adipose_Subcutaneous.h5"
# beddir <- "/media/nwknoblauch/Data/GTEx/GTEx_BED/"
# gctadir <- "/home/nwknoblauch/Downloads/gcta/gcta64"

args <- commandArgs(trailingOnly = T)

gctadir <- args[1]
h5file <- args[2]
bed_pref <-args[3]
grm_pref <- args[4]
cat_covf <- args[5]
cont_covf <- args[6]




#dgn_bed <- "/media/nwknoblauch/Data/DGN/Lev/DGN"
#gtex_h5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data.h5"
#gtex_bed <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed"
#dgn_famf <- "/media/nwknoblauch/Data/DGN/Lev/DGN.fam"
#expf <- "/media/nwknoblauch/Data/DGN/Lev/data_used_for_eqtl_study/cis_data.txt"
#gtex_famf <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed.fam"
#expdir <- "/home/nwknoblauch/Desktop/eQTL/Heritability/expdir/"
htiss <- gsub(".+/([^/]+).h5","\\1",h5file)


# bed_pref <- paste0(beddir,htiss)
bed_famf <- paste0(bed_pref,".fam")
# gds_d <- snpgdsOpen(gdsfile)
# all_sample_ids <- data_frame(whole_id=read.gdsn(index.gdsn(gds_d,"sample.id")))
# snpgdsClose(gds_d)
# snpgdsGDS2BED(gdsfile,bed_pref)

nexp <- get_rownum_h5(h5file,"EXPdata","expression")
ncov <- get_rownum_h5(h5file,"Covardat","covariates")
expmat <- read_fmat_h5(hap_h5file = h5file,groupname = "EXPdata",dataname = "expression",offset = 0,chunksize = nexp)
covmat <- read_fmat_h5(h5file,groupname="Covardat",dataname="covariates",offset=0,chunksize=ncov)
iscat <-which(apply(covmat,2,function(x){sum(round(x,1)-x)==0}))
cat_cov <- as_data_frame(covmat[,iscat])
cont_cov <- as_data_frame(covmat[,-iscat])

sexpmat <- scale(expmat,center = T,scale = T)
expleg <- read_h5_df(h5file,"EXPinfo")

expped <- read.table(bed_famf,header=F,sep="\t")
colnames(expped) <-c("Family","Individual","Paternal","Maternal","Sex",
                     "Expression")
cat_cov_df <- bind_cols(select(expped,Family,Individual),cat_cov)
cont_cov_df <- bind_cols(select(expped,Family,Individual),cont_cov)


write.table(cat_cov_df,file=cat_covf,col.names = F,row.names=F,sep="\t",quote=F)
write.table(cont_cov_df,file=cont_covf,col.names = F,row.names=F,sep="\t",quote=F)


grmc <- paste0(gctadir," --bfile ",bed_pref,
               " --maf 0.01 --make-grm --out ",
               grm_pref,"  --thread-num 10")


system(grmc)




