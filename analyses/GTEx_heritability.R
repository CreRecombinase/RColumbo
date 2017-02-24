library(gdsfmt)
library(SNPRelate)
library(RColumbo)
library(dplyr)
gtex_gdsnf <- "~/Desktop/eQTL/Whole_Blood_v6p.gds"
gtex_h5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data.h5"
gtex_bed <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed"
gtex_famf <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed.fam"
expdir <- "/home/nwknoblauch/Desktop/eQTL/Heritability/expdir/"
gtexgds <- snpgdsOpen(gtex_gdsnf)

snpgdsGDS2BED(gtexgds,gtex_bed)
nexp <- get_rownum_h5(gtex_h5,"EXPdata","orthoexpression")
expmat <- read_fmat_h5(hap_h5file = gtex_h5,groupname = "EXPdata",dataname = "orthoexpression",offset = 0,chunksize = nexp)
sexpmat <- scale(expmat,center = T,scale = T)
oexpmat <-read_fmat_h5(hap_h5file = gtex_h5,groupname = "EXPdata",dataname = "expression",offset = 0,chunksize = nexp)
expleg <- read_h5_df(gtex_h5,"EXPinfo")

expped <- read.table(gtex_famf,header=F,sep="\t")
colnames(expped) <-c("Family","Individual","Paternal","Maternal","Sex",
                     "Expression")
i <- 1


setwd("/media/nwknoblauch/Data/GTEx/")


fdfl <- list()
ndfl <- list()
expf <- paste0(expdir,"fgid_",expleg$fgeneid[i],".pheno")
for(i in i:nexp){
  cat(i," of ",nexp,"\n")

  expped <- mutate(expped,Expression=sexpmat[,i]) %>% select(Family,Individual,Expression)
  write.table(expped,file = expf,sep="\t",col.names = F,row.names = F,
              quote=F)
  gcta_command <- paste0("/home/nwknoblauch/Downloads/gcta/gcta64 --grm ",
                         "test --pheno ",expf," --reml --out ptest --thread-num 10 --reml-maxit 200")
  system(gcta_command,ignore.stdout = T,ignore.stderr = T)
  fdf <-read.table("ptest.hsq",sep="\t",nrows = 4,header=T)
  fdfl[[i]] <- fdf
  ndf <- read.table("ptest.hsq",sep="\t",skip = 5,header=F)
  ndfl[[i]] <- ndf
  file.remove(expf)
  file.remove("ptest.hsq")
}

nfdfl <- mapply(function(df,gid){
  return(mutate(df,fgeneid=gid))},fdfl,expleg$fgeneid,SIMPLIFY = F)

nndfl <- mapply(function(df,gid){
  return(mutate(df,fgeneid=gid))},ndfl,expleg$fgeneid,SIMPLIFY = F)
vpdf <- bind_rows(nfdfl)
sdf <- bind_rows(nndfl)
saveRDS(sdf,"/media/nwknoblauch/Data/GTEx/GCTA_trans_ortho_summary.RDS")
saveRDS(vpdf,"/media/nwknoblauch/Data/GTEx/GCTA_trans_ortho_estimates.RDS")

library(ggplot2)
her <- filter(vpdf,Source=="V(G)/Vp") %>% rename(h=Variance)
filter(her) %>% ggplot()+geom_histogram(aes(x=no_ortho_h),bins=100)+ggtitle("Heritability of all (orthogonalized wrt covariates) genes \n(GCTA trans+cis)")
filter(oher) %>% ggplot()+geom_histogram(aes(x=ortho_h),bins=100)+ggtitle("Heritability of all genes \n(GCTA trans+cis)")
filter(her,h>1e-06) %>% ggplot()+geom_histogram(aes(x=h),bins=100)+ggtitle("Heritability of heritable genes\n(orthogonalized wrt covariates) \n(h>1e-6) \n(GCTA trans+cis)")


ovpdf <- vpdf
osdf <- sdf

vpdf <- readRDS("/media/nwknoblauch/Data/GTEx/GCTA_trans_estimates.RDS")
sdf <- readRDS("/media/nwknoblauch/Data/GTEx/GCTA_trans_summary.RDS")

oher <- filter(ovpdf,Source=="V(G)/Vp") %>% rename(ortho_h=Variance,ortho_SE=SE) %>% select(-Source)
her <- filter(vpdf,Source=="V(G)/Vp") %>% rename(no_ortho_h=Variance,no_ortho_SE=SE) %>% select(-Source)
nher <- inner_join(oher,her,by="fgeneid")
ggplot(nher)+geom_point(aes(x=ortho_h,y=no_ortho_h))+ggtitle("Heritability of raw expression \nvs\n heritability of expression orthogonalized to covariates")





