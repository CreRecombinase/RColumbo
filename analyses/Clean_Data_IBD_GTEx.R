#Clean data for
library(RColumbo)
library(rhdf5)
library(dplyr)
library(tidyr)
rawh5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_raw_data.h5"
haph5 <- "/home/nwknoblauch/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno_hap.h5"
haph5s <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
gwash5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/IBD.h5"
dbsnph5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/dbsnp.h5"
oh5s <- paste0("/media/nwknoblauch/Data/GTEx/chr",1:22,"_1kg_Whole_Blood_IBD.h5")


gwas_df <- read_h5_df(gwash5,"GWAS") %>% mutate(gwas_ind=1:length(rsid))
dbsnp_df <- read_h5_df(dbsnph5,"dbSNP")
both_df <- inner_join(gwas_df,dbsnp_df)
gtex_df <- read_h5_df(rawh5,"SNPinfo") %>% mutate(gtex_ind=1:length(rsid))
ngtex_df <- semi_join(gtex_df,both_df) %>% filter(!duplicated(rsid))
nboth_df <- semi_join(both_df,ngtex_df) %>% filter(!duplicated(rsid))
fgtex_df <- inner_join(ngtex_df,nboth_df)
fgtex_df <- select(fgtex_df,-doFlip,-N)
#int write_Rnumeric_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::NumericVector &data,const unsigned int deflate_level);
i <- 1

for(i in 2:22){
  cat(paste0(i,"of 22\n"))
  haph5 <- haph5s[i]
  oh5 <- oh5s[i]
  tfgtex <- filter(fgtex_df,chrom==i)
  legdf <- read_h5_df(haph5,"Legend") %>% mutate(hap_ind=1:length(rsidi),chrom=i)
  mapdf <- read_h5_df(haph5,"Map") %>% mutate(map_ind=1:length(rsidi),chrom=i)
  legmapdf <- inner_join(legdf,mapdf) %>% select(-rsid) %>% rename(rsid=rsidi)
  nlegmapdf <- inner_join(legmapdf,tfgtex)
  thapdat <- flip_hap(haph5,nlegmapdf$hap_ind,0,length(nlegmapdf$hap_ind),length(nlegmapdf$hap_ind))
  tgtexdat <- read_dmat_chunk_ind(rawh5,"SNPdata","orthogenotype",nlegmapdf$gtex_ind)
  wlegmapdf <- select(nlegmapdf,chrom,pos,rsid,cummap,beta,serr)
  write_h5_df(wlegmapdf,"Legend",oh5)
  write_dmatrix_h5(h5file = oh5,groupname="Haplotype",dataname="genotype",Nsnps = ncol(thapdat),Nind = nrow(thapdat),data = thapdat,deflate_level = 2)
  write_dmatrix_h5(h5file = oh5,groupname="eQTL",dataname="genotype",Nsnps = ncol(tgtexdat),Nind = nrow(tgtexdat),data = tgtexdat,deflate_level = 2)
  H5close()
  # write_Rnumeric_h5(h5file = oh5,groupname = "Haplotype",dataname = "genotype",thapdat,deflate_level = 2)
  # write_Rnumeric_h5(h5file = oh5,groupname = "eQTL",dataname = "genotype",tgtexdat,deflate_level = 2)
  gc()
}

