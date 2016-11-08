library(RColumbo)
library(dplyr)
library(tidyr)
library(readr)
library(rhdf5)
library(gdsfmt)


setwd("/media/nwknoblauch/Data/DGN/Lev/")
dgn_covarbf <- "covariates/Biological_and_hidden_factors.txt"
dgn_covartf <- "covariates/Technical_factors.txt"
f_DGN_geno_h5 <- "eQTL/DGN_geno_exp.h5"
ortho_DGN <- "eQTL/DGN_ortho.h5"
dgn_expf <- "data_used_for_eqtl_study/cis_iso_data.txt"
dgn_gdsf <- "DGN.gds"



dgn_gds <- gdsfmt::openfn.gds(dgn_gdsf)



iso_dat <- read.table(dgn_expf,sep="\t",header=T,stringsAsFactors = F)
iso_dat <- rename(iso_dat,Id=X)
iso_id <- as.character(iso_dat$Id)
iso_dat <- data.matrix(select(iso_dat,-Id))
iso_col <- data_frame(isocol=colnames(iso_dat))
iso_col <- mutate(iso_col,coord=gsub(".+(chr.+)","\\1",isocol))
niso_col <- separate(iso_col,coord,into = c("chr","start","stop"))
niso_col <-mutate(niso_col,chr=as.integer(gsub("chr([0-9]+)","\\1",chr)))
niso_col <- mutate(niso_col,fgeneid=1:n())

niso_col <- filter(niso_col,!is.na(chr))

exp_leg <- select(niso_col,-isocol) %>% mutate(start=as.integer(start),stop=as.integer(stop),chr=as.integer(chr)) %>% rename(end=stop) %>% distinct(chr,start,end,.keep_all=T)
iso_dat <- iso_dat[,exp_leg$fgeneid]



dgn_leg <- data_frame(ind=read.gdsn(index.gdsn(dgn_gds,"snp.id")),                     chrom=as.integer(read.gdsn(index.gdsn(dgn_gds,"snp.chromosome"))),
                      pos=as.integer(read.gdsn(index.gdsn(dgn_gds,"snp.position"))))

snp_leg <- select(dgn_leg,chrom,pos) %>% mutate(ind=1:n()) %>% distinct(chrom,pos,.keep_all=T)



dgn_covarb <- read.table(dgn_covarbf,header=T,sep="\t",stringsAsFactors = F)
dgn_covart <- read.table(dgn_covartf,header=T,sep="\t",stringsAsFactors = F)


dgn_snp <- read.gdsn(index.gdsn(dgn_gds,"genotype"))
dgn_ids <- read.gdsn(index.gdsn(dgn_gds,"sample.id"))
dgn_ids <-gsub("WG[0-9]+-DNA.+_[0-9]+_(LD[0-9]+)-.+","\\1",x = dgn_ids)

rownames(dgn_snp) <- dgn_ids
dgn_snp <- dgn_snp[iso_id,]
dgn_snp <- dgn_snp[,snp_leg$ind]
total_covar <- inner_join(dgn_covarb,dgn_covart,by="X")
total_covmat <-data.matrix(select(total_covar,-X))

stopifnot(nrow(iso_dat)==nrow(total_covar),nrow(iso_dat)==nrow(dgn_snp))


if(file.exists(f_DGN_geno_h5)){
  file.remove(f_DGN_geno_h5)
}

write_dmatrix_h5(f_DGN_geno_h5,"SNPdata","genotype",
                 Nsnps = ncol(dgn_snp),
                 Nind = nrow(dgn_snp),
                 data = dgn_snp,deflate_level = 4)

write_dmatrix_h5(f_DGN_geno_h5,"EXPdata","expression",
                 Nsnps=ncol(iso_dat),
                 Nind=nrow(iso_dat),
                 data=data.matrix(iso_dat),deflate_level=4)


write_dmatrix_h5(f_DGN_geno_h5,"Covardat","covariates",
                 Nsnps = ncol(total_covmat),Nind = nrow(total_covmat),data = total_covmat,deflate_level = 4)

if(file.exists(ortho_DGN)){
  file.remove(ortho_DGN)
}

orthogonalize_dataset(h5filename = f_DGN_geno_h5,newh5filename = ortho_DGN,
                      datagroup = "SNPdata",datasetname = "genotype",chunksize = 30000,deflate_level = 4,
                      newdatasetname = "genotype")
orthogonalize_dataset(h5filename = f_DGN_geno_h5,newh5filename = ortho_DGN,
                      datagroup = "EXPdata",datasetname = "expression",chunksize = 3000,deflate_level = 4,
                      newdatasetname = "expression")
write_h5_df(snp_leg,group = "SNPinfo",outfile = ortho_DGN,deflate_level = 4)
write_h5_df(df=exp_leg,group = "EXPinfo",outfile = ortho_DGN,deflate_level=4)



chroms <- 2:22
for(i in chroms){
  chromosome <- i
  outh5file <- paste0("eQTL/DGN_Chr",i,"_ortho_noflip.h5")
  run_eqtl(rawh5 = ortho_DGN,outh5file,chromosome,chunksize=30000,cis_pcutoff=0.05,trans_pcutoff=5e-5,cisdist_cutoff=1e6,append=F,useortho=F)
}
