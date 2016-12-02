#Code to prep for Lambda Estimation
library(readr)
library(ggplot2)
library(dplyr)
library(gdsfmt)
library(SNPRelate)
library(RColumbo)
library(feather)

# gtexf <- "~/Desktop/eQTL/Snake/WholeBlood.txt"
gtex_gdsnf <- "~/Desktop/eQTL/Whole_Blood_v6p.gds"



gtex_h5 <- "~/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data_flip_ortho.h5"
geno_h5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data_flip.h5"
idgtexf <- "/media/nwknoblauch/Data/GTEx/ind_gtex_snps.RDS"
ind_gtex_gwasf <- "/media/nwknoblauch/Data/GTEx/ind_gtex_gwas.RDS"
eqtlh5s <- paste0("/media/nwknoblauch/Data/GTEx/GTEx_v6p_h5files/ortho_flip/WB_Chr",1:22,"_v6p_ortho_flip.h5")
cisresf <- "/media/nwknoblauch/Data/GTEx/cis_res_ortho_flip.RDS"
transresf <- "/media/nwknoblauch/Data/GTEx/trans_res_ortho_flip.RDS"
transresfeather <- "/media/nwknoblauch/Data/GTEx/trans_res_ortho_flip.feather"
allresfeather <- "/media/nwknoblauch/Data/GTEx/all_res_ortho_flip.feather"

eqtl_legf <- "/home/nwknoblauch/Desktop/eQTL/Snake/eqtl_leg.RDS"
gwas_f <- "/media/nwknoblauch/Data/gwas_data/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc.gz"
gwasRDS <- "/media/nwknoblauch/Data/gwas_data/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc.RDS"

overwrite=TRUE
# install_github("mannau/h5",args="--configure-vars='LIBS=-L/usr/local/lib  -lhdf5 -lhdf5_hl -lhdf5_cpp -lz -lblosc CPPFLAGS=-I/usr/local/include '")

gtex_leg <- read_h5_df(geno_h5,"SNPinfo") %>% mutate(snp_ind=1:n())
group_by(gtex_leg,chrom) %>% summarise(snpct=n()) %>% arrange(desc(snpct)) %>% data.frame()
saveRDS(gtex_leg,eqtl_legf)

chleg <- filter(gtex_leg,chrom==19)


teqtlh5 <- eqtlh5s[2]
bcdf <- read_h5_df(teqtlh5,groupname = "cis_eQTL",subcols = c("theta","chrom"))

library(RColumbo)
library(dplyr)
eqtlh5s <- paste0("WB_Chr",1:22,"_v6p_ortho_flip.h5")
# gwasRDS <- "EUR.IBD.gwas.assoc.RDS"
# cisresf <- "cis_res_ortho_flip.RDS"
# transresf <- "trans_res_ortho_flip.RDS"
gwas_dat <- readRDS(gwasRDS)
cis_res <- bind_rows(lapply(eqtlh5s,sample_h5_df,blocksize=20000,pval.cutoff=0.05,groupname="cis_eQTL",gwasdf=gwas_dat))
cis_res <- mutate(cis_res,isCis=T)
 saveRDS(cis_res,cisresf)

 trans_res <- bind_rows(lapply(eqtlh5s,sample_h5_df,blocksize=20000,pval.cutoff=0.005,groupname="trans_eQTL",gwasdf=gwas_dat))
trans_res <- mutate(trans_res,isCis=F)
system.time(saveRDS(trans_res,transresf))
system.time(write_feather(trans_res,transresfeather))


###Import and transform GWAS data (then save to RDS)
gwasf2 <- "/media/nwknoblauch/Data/gwas_data/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc_2.feather"
if(!file.exists(gwasRDS)|overwrite){
  gwas_dat <- read_delim(gwas_f,delim="\t",col_names=T)
  gwas_dat <- rename(gwas_dat,chrom=CHR,rsid=SNP,pos=BP,ref_allele=A2,eff_allele=A1,p.value=P,serr.g=SE)
  gwas_dat <- mutate(gwas_dat,posl=strsplit(Direction,split=""))%>% mutate(npos=map_int(posl,~sum(.=="+")),
                                                                           nneg=map_int(posl,~sum(.=="-")),
                                                                           namb=map_int(posl,~sum(.=="?"))) %>% dplyr::select(-posl)
  gwas_dat <- filter(gwas_dat,npos+nneg==15)
  gwas_dat <- select(gwas_dat,chrom,rsid,pos,eff_allele,ref_allele,OR,serr.g,p.value)
  write_feather(gwas_dat,gwasf2)
  gwas_dat <- mutate(gwas_dat,doFlip=as.integer(A1>A2),
                     nOR=ifelse(doFlip==1,1/OR,OR),posl=strsplit(Direction,split=""))%>% mutate(npos=map_int(posl,~sum(.=="+")),
                                                                                                nneg=map_int(posl,~sum(.=="-")),
                                                                                                namb=map_int(posl,~sum(.=="?"))) %>% dplyr::select(-posl)

  gwas_dat <- filter(gwas_dat,npos+nneg==15)
  gwas_dat <- select(gwas_dat,chrom=CHR,rsid=SNP,pos=BP,beta=nOR,serrg=SE,pvalg=P) %>% mutate(beta=log(beta))
  saveRDS(gwas_dat,gwasRDS)
}else{
  gwas_dat <- readRDS(gwasRDS)
}


#
# ###Import and transform GTEx Genotype data (then save independent SNPs to RDS)
# if((!file.exists(gtex_gdsnf)|(!file.exists(idgtexf)))){
#   if(!file.exists(gtex_gdsnf)){
#   worked <- write_gtex_gdsn(gtexf,gtex_gdsnf,chunksize = 500000)
#   }
#   gtexgds <- snpgdsOpen(gtex_gdsnf)
#
#
#
#   idgtex <- snpgdsLDpruning(gtexgds,autosome.only=T)
#   idgtexdf <- bind_rows(mapply(function(ch,id){
#     data_frame(chrom=ch,ind=id)
#   },1:22,idgtex,SIMPLIFY = F))
#   gtex_leg <- data_frame(ind=as.integer(read.gdsn(index.gdsn(gtexgds,"snp.id"))),
#                          chrom=as.integer(read.gdsn(index.gdsn(gtexgds,"snp.chromosome"))),
#                          pos=as.integer(read.gdsn(index.gdsn(gtexgds,"snp.position"))))
#
#
#   ind_gtex_leg <- semi_join(gtex_leg,idgtexdf) %>% select(-ind)
#   saveRDS(ind_gtex_leg,idgtexf)
#   rm(gtex_leg,idgtex,idgtexdf,worked)
# }else{
#   ind_gtex_leg <- readRDS(idgtexf)
# }


##We now also need a method for getting independent eQTL that is not agnostic to eQTL effect size.
## This is especially important for trans, where signal is sparse)
##Let's see how fast it is to do one gene and one chromosome at a time (maybe we can just throw it on the cluster)





# if(!file.exists(ind_gtex_gwasf)){
#   ind_gwas_dat <- inner_join(gwas_dat,ind_gtex_leg,by=c("chrom","pos")) %>% dplyr::select(-rsid)
#   saveRDS(ind_gwas_dat,ind_gtex_gwasf)
# }else{
#   ind_gwas_dat <- readRDS(ind_gtex_gwasf)
# }


if(!file.exists(cisresf)|overwrite){
  cis_res <- bind_rows(lapply(eqtlh5s,fsample_h5_df,blocksize=20000,pval.cutoff=0.05,groupname="cis_eQTL",gwasdf=gwas_dat))
  cis_res <- mutate(cis_res,isCis=T)
  saveRDS(cis_res,cisresf)
  # write.table(cis_res,file = "/media/nwknoblauch/Data/GTEx/cis_res_no_ortho.txt",sep="\t",col.names=T,row.names=F,quote=F)
}else{
  cis_res <- readRDS(cisresf)
}

if(!file.exists(transresf)|overwrite){
  trans_res <- bind_rows(lapply(eqtlh5s,fsample_h5_df,blocksize=20000,pval.cutoff=0.005,groupname="trans_eQTL",gwasdf=gwas_dat))
  trans_res <- mutate(trans_res,isCis=F,tstat=theta/serr)
  saveRDS(trans_res,transresf)
}else{
  trans_res <- readRDS(transresf)
}

all_res <- bind_rows(cis_res,trans_res)

if(!file.exists(allresfeather)){
  all_res <- bind_rows(readRDS(transresf),readRDS(cisresf))
  # saveRDS(all_res,allresf)
  write_feather(all_res,allresfeather)
}

ggplot(all_res)+geom_histogram(aes(x=beta))









