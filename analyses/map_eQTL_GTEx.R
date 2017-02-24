library(RColumbo)
library(testthat)
library(dplyr)
library(tidyr)
library(readr)
library(rhdf5)


wbsnpf <- "~/Desktop/eQTL/Snake/WholeBlood.txt"
wbexpf <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_expression_matrices/Whole_Blood_Analysis.v6p.normalized.expression.bed.gz"

covarf <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_covariates/Whole_Blood_Analysis.v6p.covariates.txt"
ibdh5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/IBD.h5"
gwasRDS <- "/media/nwknoblauch/Data/gwas_data/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc.RDS"
h5file <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data_no_flip.h5"
oh5file <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data_noflip_ortho.h5"
exp_legdf <- read_h5_df(oh5file,"EXPinfo") %>% arrange(chrom,start,end)
snp_legdf <- read_h5_df(oh5file,"SNPinfo") %>% arrange(chrom,pos)

snp_ldfl <- split(snp_legdf,snp_legdf$chrom)
exp_ldfl <- split(exp_legdf,exp_legdf$chrom)






tsnp <- snp_ldfl[[1]]
texp <- exp_ldfl[[1]]

ttsnp <- chunk_df(tsnp,chunk.size = 10000)


if(file.exists(h5file)){
file.remove(h5file)
}
if(file.exists(oh5file)){
  file.remove(oh5file)
}
read_gtex_snp(wbsnpf,h5file,chunksize=5e5,FlipAllele = F)
read_gtex_expression(gzfile = wbexpf,h5file = h5file)
write_covar_h5(covarf=covarf,h5file = h5file,chunksize = 1,deflate_level = 4)

orthogonalize_dataset(h5file,oh5file,"SNPdata","genotype","genotype",500000,4)
orthogonalize_dataset(h5file,oh5file,"EXPdata","expression","expression",5000,4)
snpleg <- read_h5_df(h5file,"SNPinfo")
write_h5_df(snpleg,"SNPinfo",oh5file)
expleg <- read_h5_df(h5file,"EXPinfo")
write_h5_df(expleg,"EXPinfo",oh5file)

snpdata <- read_fmat_h5(oh5file,groupname = "SNPdata",dataname="genotype",offset=0,chunksize=get_rownum_h5(oh5file,groupname=""))


#dbsnpdf <- read_h5_df(dbsnph5,"dbSNP")
eqtlleg <- read_h5_df(oh5file,"SNPinfo") %>% mutate(eqtl_ind=1:n())
ibdleg <- readRDS(gwasRDS)

eqtl_rsid <- inner_join(ibdleg,eqtlleg,by=c("chrom","pos"))
group_by(eqtlleg,chrom) %>% summarise(nsnp=n()) %>% arrange(desc(nsnp))


chroms <- 1:22
oh5file <- "Whole_Blood_eQTL_v6p_raw_data_flip_ortho.h5"



oh5file <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data_flip_ortho.h5"
for(i in chroms){
  chromosome <- i
  outh5file <- paste0("/media/nwknoblauch/Data/GTEx/GTEx_v6p_h5files/WB_Chr",i,"_v6p_ortho_flip.h5")
  run_eqtl(rawh5 = oh5file,outh5file,chromosome,chunksize=30000,cis_pcutoff=0.05,trans_pcutoff=5e-5,cisdist_cutoff=1e6,append=F,useortho=F)
}




data_h5files <-paste0("/home/nwknoblauch/Desktop/eQTL/GTEx/GTEx_WB_Chr",1:22,"_no_flip_orth.h5")
snpleg <- read_h5_df(h5file,"SNPinfo")
snpleg <- mutate(snpleg,snp_ind=1:n())
expleg <- read_h5_df(h5file,"EXPinfo") %>% mutate(exp_ind=1:n()) %>% rename(chrom=chr)

chroms <- 1:22
for(i in chroms){
  cat(i)
  tsnpleg <- filter(snpleg,chrom==i)
  texpleg <- filter(expleg,chrom==i)
  tsnp <- read_fmat_chunk_ind(h5file = oh5file,groupname = "SNPdata","genotype",indvec = tsnpleg$snp_ind)
  write_h5_df(data_h5files[i],group ="SNPinfo",df = tsnpleg,deflate_level = 4)
  write_dmatrix_h5(data_h5files[i],groupname = "SNPdata",dataname = "genotype",Nsnps = ncol(tsnp),Nind = nrow(tsnp),data = tsnp,deflate_level = 4)
  rm(tsnp)
  gc()
  write_h5_df(data_h5files[i],group="EXPinfo",df = texpleg,deflate_level = 4)
  texp <- read_fmat_chunk_ind(h5file = oh5file,group = "EXPdata","expression",indvec = texpleg$exp_ind)
  write_dmatrix_h5(data_h5files[i],group = "EXPdata",dataname = "expression",Nsnps = ncol(texp),Nind = nrow(texp),data = texp,deflate_level = 4)
  rm(texp)
  gc()
}

