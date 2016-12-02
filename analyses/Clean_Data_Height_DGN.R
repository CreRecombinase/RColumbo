#Clean data for
library(RColumbo)
library(rhdf5)
library(dplyr)
library(tidyr)


gwdh <- strsplit(getwd(),split = "/",fixed=T)[[1]][3]
if(gwdh=="nwknoblauch"){
  rawh5 <- "/media/nwknoblauch/Data/misc_h5/DGN.h5"
  haph5s <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
  oh5s <- paste0("/media/nwknoblauch/Data/DGN/EUR.chr",1:22,"_1kg_DGN_Height.h5")
  gwasgz <- "~/Desktop/eQTL/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"
  # oh5files <- paste0("/media/nwknoblauch/Data/GTEx/LDmats/chr",1:22,"_1kg_WB_Height_LD.h5")
}else{

  h5files <- paste0("/group/xhe-lab/RSS/EUR.chr",1:22,"_1kg_Whole_Blood_HEIGHT.h5")
  oh5files <- paste0("/group/xhe-lab/RSS/LD_matrices/chr",1:22,"_1kg_WB_Height_LD.h5")
}
raw_legend_group <- "Legend"
raw_data_group <- "Genotype"
raw_data_data <- "genotype"
# rawh5 <- "/media/nwknoblauch/Data/misc_h5/DGN.h5"
# haph5 <- "/media/nwknoblauch/Data/GTEx/chr22_1kg_Whole_Blood_HEIGHT.h5"


# haph5s <- paste0("/group/xhe-lab/RSS/IMPUTE/EUR.chr",1:22,"_1kg_geno_hap.h5")
# gwasgz <-
# gwasgz <- "/group/xhe-lab/RSS/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"

#gwash5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/IBD.h5"
#dbsnph5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/dbsnp.h5"



#gwas_df <- read_h5_df(gwash5,"GWAS") %>% mutate(gwas_ind=1:length(rsid))
gwas_df <- read.table(gwasgz,header=T,stringsAsFactors = F) %>% mutate(beta=ifelse(Allele1<Allele2,-b,b))
gwas_df <- rename(gwas_df,rsid=MarkerName) %>% mutate(gwas_ind=1:n())
                                        #dbsnp_df <- read_h5_df(dbsnph5,"dbSNP") %>% mutate(rsid=paste0("rs",rsid))

                                        #both_df <- inner_join(gwas_df,dbsnp_df)
both_df <- gwas_df
both_df <- mutate(both_df,rsid=as.integer(gsub("rs","",rsid)))
both_df <- filter(both_df,!is.na(rsid))
gtex_df <- read_h5_df(rawh5,raw_legend_group) %>% mutate(raw_ind=1:length(rsid))

ngtex_df <- semi_join(gtex_df,both_df) %>% filter(!duplicated(rsid),!is.na(rsid))
nboth_df <- semi_join(both_df,ngtex_df) %>% filter(!duplicated(rsid),!is.na(rsid))
fgtex_df <- inner_join(ngtex_df,nboth_df)
fgtex_df <- select(fgtex_df,-doFlip,-N)
#int write_Rnumeric_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::NumericVector &data,const unsigned int deflate_level);
i <- 22

for(i in 1:21){
  cat(paste0(i,"of 22\n"))
  haph5 <- haph5s[i]
  oh5 <- oh5s[i]
  tfgtex <- filter(fgtex_df,chrom==i)
  legdf <- read_h5_df(haph5,"Legend") %>% mutate(hap_ind=1:n(),chrom=i) %>% mutate(rsidi=as.integer(gsub("rs","",rsid)))
  mapdf <- read_h5_df(haph5,"Map") %>% mutate(map_ind=1:n(),chrom=i)
  legmapdf <- inner_join(legdf,mapdf) %>% select(-rsid) %>% rename(rsid=rsidi)
  nlegmapdf <- inner_join(legmapdf,tfgtex)
  thapdat <- flip_hap(haph5,nlegmapdf$hap_ind,0,length(nlegmapdf$hap_ind),length(nlegmapdf$hap_ind))
  tgtexdat <- read_dmat_chunk_ind(rawh5,raw_data_group,raw_data_data,nlegmapdf$raw_ind)
  wlegmapdf <- select(nlegmapdf,chrom,pos,rsid,cummap,beta,serr=SE)
  write_h5_df(wlegmapdf,"Legend",oh5)
  write_dmatrix_h5(h5file = oh5,groupname="Haplotype",dataname="genotype",Nsnps = ncol(thapdat),Nind = nrow(thapdat),data = thapdat,deflate_level = 2)
  write_dmatrix_h5(h5file = oh5,groupname="eQTL",dataname="genotype",Nsnps = ncol(tgtexdat),Nind = nrow(tgtexdat),data = tgtexdat,deflate_level = 2)
  H5close()
  # write_Rnumeric_h5(h5file = oh5,groupname = "Haplotype",dataname = "genotype",thapdat,deflate_level = 2)
  # write_Rnumeric_h5(h5file = oh5,groupname = "eQTL",dataname = "genotype",tgtexdat,deflate_level = 2)
  gc()
}

