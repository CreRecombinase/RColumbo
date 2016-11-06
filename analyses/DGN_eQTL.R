library(RColumbo)
library(dplyr)
library(tidyr)
library(readr)
library(rhdf5)
library(gdsfmt)


setwd("/project/xinhe/Lev")
dgn_covarbf <- "covariates/Biological_and_hidden_factors.txt"
dgn_covartf <- "covariates/Technical_factors.txt"
f_DGN_geno_h5 <- "eQTL/DGN_geno_exp.h5"
ortho_DGN <- "eQTL/DGN_ortho.h5"
dgn_expf <- "data_used_for_eqtl_study/cis_iso_data.txt"
dgn_gdsf <- "DGN.gds"



dgn_gds <- gdsfmt::openfn.gds(dgn_gdsf)



iso_dat <- read.table(dgn_expf,sep="\t",header=T)
iso_dat <- rename(iso_dat,Id=X)
iso_col <- data_frame(isocol=colnames(iso_dat)[-1])
iso_col <- mutate(iso_col,coord=gsub(".+(chr.+)","\\1",isocol))
niso_col <- separate(iso_col,coord,into = c("chr","start","stop"))
niso_col <-mutate(niso_col,chr=as.integer(gsub("chr([0-9]+)","\\1",chr)))
niso_col <- mutate(niso_col,fgeneid=1:n())
iso_dat <- iso_dat[,c(T,!is.na(niso_col$chr))]
niso_col <- filter(niso_col,!is.na(chr))

exp_leg <- select(niso_col,-isocol) %>% mutate(start=as.integer(start),stop=as.integer(stop),chr=as.integer(chr)) %>% rename(end=stop)
niso_col <- filter(niso_col,!is.na(chr))



dgn_leg <- data_frame(ind=read.gdsn(index.gdsn(dgn_gds,"snp.id")),                     chrom=as.integer(read.gdsn(index.gdsn(dgn_gds,"snp.chromosome"))),
                      pos=as.integer(read.gdsn(index.gdsn(dgn_gds,"snp.position"))))

dgn_leg <- select(dgn_leg,chrom,pos) %>% mutate(ind=1:n())



write_h5_df(df=exp_leg,group = "EXPinfo",outfile = ortho_DGN)
write_h5_df(df=dgn_leg,group = "SNP")

n_iso_dat <- dplyr::select(iso_dat,one_of(c("Id",niso_col$isocol)))

n_iso_dat <-select(n_iso_dat,-Id)
n_iso_dat <- n_iso_dat[,niso_col$isocol]





anti_join(dgn_col,select(iso_col,gene=gene2)) %>% dim()


iso_col <- separate(iso_col,isocol,into=c("NM","entrez","symbol","symbol","TSS","loc"),sep="_",remove = F)


dgn_exp <- read.table(dgn_expf,header=T,sep = "\t")
dgn_exp <- rename(dgn_exp,Id=X) %>% select(-X.1)

dgn_col <- data_frame(gene=colnames(dgn_exp)[-1])
dgn_covarb <- read.table(dgn_covarbf,header=T,sep="\t")
dgn_covart <- read.table(dgn_covartf,header=T,sep="\t")






dgn_snp <- read.gdsn(index.gdsn(dgn_gds,"genotype"))
dgn_ids <- read.gdsn(index.gdsn(dgn_gds,"sample.id"))







t_dgn_ids <- gsub("WG[0-9]+-DNA.+_[0-9]+_(LD[0-9]+)-.+","\\1",x = dgn_ids)
rownames(dgn_snp) <- t_dgn_ids
dgn_snp <- dgn_snp[dgn_exp$Id,]
total_covar <- inner_join(dgn_covarb,dgn_covart,by="X")
total_covmat <-data.matrix(select(total_covar,-X))
dgn_expmat <- data.matrix(select(dgn_exp,-Id))


write_dmatrix_h5(f_DGN_geno_h5,"SNPdata","genotype",Nsnps = ncol(dgn_snp),Nind = nrow(dgn_snp),data = dgn_snp,deflate_level = 4)

write_dmatrix_h5(f_DGN_geno_h5,"EXPdata","expression",Nsnps=ncol(dgn_expmat),Nind=nrow(dgn_expmat),data=data.matrix(dgn_expmat),deflate_level=4)


write_dmatrix_h5(f_DGN_geno_h5,"Covardat","covariates",
                 Nsnps = ncol(total_covmat),Nind = nrow(total_covmat),data = total_covmat,deflate_level = 4)

orthogonalize_dataset(h5filename = f_DGN_geno_h5,newh5filename = ortho_DGN,
                      datagroup = "SNPdata",datasetname = "genotype",chunksize = 30000,deflate_level = 4,
                      newdatasetname = "genotype")
orthogonalize_dataset(h5filename = f_DGN_geno_h5,newh5filename = ortho_DGN,
                      datagroup = "EXPdata",datasetname = "expression",chunksize = 3000,deflate_level = 4,
                      newdatasetname = "expression")
dgn_leg <-select(dgn_leg,chrom,pos)
write_h5_df(dgn_leg,group = "SNPinfo",outfile = ortho_DGN,deflate_level = 4)
write_h5_df(df=exp_leg,group = "EXPinfo",outfile = ortho_DGN,deflate_level=4)



chroms <- 2:22
oh5file <- "Whole_Blood_eQTL_v6p_raw_data_flip_ortho.h5"

for(i in chroms){
  chromosome <- i
  outh5file <- paste0("/media/nwknoblauch/Data/DGN/Lev/eQTL/DGN_Chr",i,"_ortho_noflip.h5")
  run_eqtl(rawh5 = ortho_DGN,outh5file,chromosome,chunksize=30000,cis_pcutoff=0.05,trans_pcutoff=5e-5,cisdist_cutoff=1e6,append=F,useortho=F)
}
