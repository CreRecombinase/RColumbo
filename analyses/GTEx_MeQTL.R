#Code to try to reconstruct eQTL
library(RColumbo)
library(rhdf5)
library(dplyr)
library(tidyr)
eqtlh5 <- "~/Desktop/eQTL/Snake/IBD_WholeBlood_eQTL.h5"
rawh5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_raw_data.h5"
haph5 <- "/home/nwknoblauch/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno_hap.h5"
haph5s <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
gwash5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/IBD.h5"
oh5s <- paste0("/media/nwknoblauch/Data/GTEx/chr",1:22,"_Whole_Blood_eQTL.h5")

gwas_eqtl_snps <- intersect_col(h5file1 = gwash5,"GWAS","rsid",rawh5,"SNPinfo","rsid")
for(i in 2:22){
  cat(paste0(i," of 22\n"))
  rewrite_h5(haph5s[i])
}


haplotype_snps <- list()
for(i in 2:22){
  cat(paste0(i,"of 22\n"))
  haplotype_snps[[i]] <- intersect(c(h5read(haph5s[i],"/Legend/rsid")),c(h5read(haph5s[i],"/Map/rsid")))
}


snpleg <- read_h5_df(rawh5,"SNPinfo") %>% mutate(snp_ind=1:length(rsid))
expleg <- as_data_frame(t(h5read(rawh5,"EXPinfo/annomat"))) %>% rename(fgeneid=V1,sgeneid=V2,chrom=V3,TSStart=V4,TSStop=V5) %>% mutate(exp_id=1:(length(fgeneid)))
#sub_snpleg <- subset_all_hap(snpleg,haph5s,gwash5)
#saveRDS(sub_snpleg,"~/Desktop/eQTL/hap_eQTL.RDS")
sub_snpleg <- readRDS("~/Desktop/eQTL/hap_eQTL.RDS")
rsidvec <- as.integer(gsub("rs","",sub_snpleg$rsid))




for(i in 2:22){
  cat(paste0(i," of ",22,"\n"))
  oh5 <- oh5s[i]
  haph5 <- haph5s[i]
  stat_extract(eqtlh5,rawh5,haph5,gwash5,oh5,i,15000,cis_pcutoff=0.01,
               trans_pcutoff=1e-3,cis_LDcutoff=0.8,trans_LDcutoff=0.8,cisdist_cutoff=1e6,append=F)
}

gwasd <- read_h5_df(gwash5,"GWAS")

siggwasd <- filter(gwasd,-log(pval,base=10)>8) %>% arrange(pval) %>% select(-N)

g_pos <- mutate(gwasd,rsid=paste0("rs",rsid)) %>% inner_join(sub_snpleg,by="rsid")
sort(table(sig_pos$chrom))


cis_eqtl1 <- read_h5_df(oh5s[1],"cis_eQTL")
cis_eqtl1 <- mutate(cis_eqtl1,ensid=factor(paste0("ENSG",sprintf("%011d",fgeneid)))) %>% rename(see=serr)
trans_eqtl1 <- read_h5_df(oh5s[1],"trans_eQTL")
gwas_cis <- mutate(cis_eqtl1,rsid=paste0("rs",rsid)) %>%inner_join(g_pos,by=c("rsid"))
#snpdf <- read_h5_df(eqtlh5,"SNP")
eqtldf <- read_h5_df(eqtlh5,"eQTL")



gc()
gmd <-vegas.merge(vegasf,gmf)
gmd <- filter(gmd,!duplicated(ensid))
gmd <- mutate(gmd,fgeneid=as.integer(gsub("ENSG","",ensid)))

eqtl_vegas <- inner_join(gwas_cis,gmd,by="ensid") %>% arrange(vegp) %>% select(-cummap,-snp_ind,-allele0,-allele1,-doFlip,-hsnpid,-N) %>% mutate(isCis=1)

ilr_cis <- filter(eqtl_vegas,ensid=="ENSG00000162594")
ilrfgid <- ilr_cis$fgeneid[1]
acis <-filter(eqtldf,fgeneid==ilrfgid)
e_cis_rss <- select(acis,fgeneid,rsid,thetahat,see) %>% mutate(isCis=1)


ilr_trans <- filter(trans_eqtl1,fgeneid==ilrfgid) %>% mutate(isCis=0)
e_trans_rss <- select(ilr_trans,fgeneid,rsid,thetahat,see=serr,isCis)
e_rss <- bind_rows(e_cis_rss,e_trans_rss)
df <- select(e_rss,rsid,betahat=thetahat,serr=see) %>% mutate(rsid=paste0("rs",rsid),Nind=338)
df <- semi_join(df,sub_snpleg,by="rsid")
head(sub_snpleg)



ilr_trans <- mutate(ilr_trans,rsid=paste0("rs",rsid)) %>% rename(see=serr) %>% inner_join(g_pos,by="rsid")
ilr_trans <- select(ilr_trans,-N,-allele0,-allele1,-doFlip,-hsnpid,-cummap,-snp_ind)
ilr_filter <- function(x){x==ilrfgid}
trans_dfl <- list()
for(i in 1:22){
  cat(paste0(i," of 22\n"))
  trans_dfl[[i]] <- read_h5_df_filter(oh5s[i],groupname = "trans_eQTL",colname = "fgeneid",filter_fun = ilr_filter)
}
nilr <- bind_rows()


all_ilr <- bind_rows(trans_dfl)
all_ilr <- mutate(all_ilr,isCis=0)
all_ilr <- mutate(all_ilr,rsid=paste0("rs",rsid)) %>% rename(see=serr) %>% inner_join(g_pos,by="rsid")
all_ilr <- select(all_ilr,-N,-allele0,-allele1,-doFlip,-hsnpid,-cummap,-snp_ind)
all_ilr <- mutate(all_ilr,ensid=factor(paste0("ENSG",sprintf("%011d",fgeneid))))
trans_ilr <-inner_join(all_ilr,gmd,by="ensid") %>% arrange(vegp)

vegas_ilr <- bind_rows(trans_ilr,ilr_cis)
vegas_ilr <- mutate(vegas_ilr,eqtlt=thetahat/see)
group_by(vegas_ilr,chrom) %>% summarise(nsig=sum(abs(eqtlt)<3.5),tot=n(),meannsig=nsig/tot) %>% arrange(meannsig)


preddbf <- "~/Downloads/PredictDB_Covariance/TW_Whole_Blood/TW_Whole_Blood_0.5.db"
preddb <- src_sqlite(preddbf,create=F)
predw <- tbl(preddb,"weights")
predw <- collect(predw,n=Inf) %>% mutate(ensid=gsub("\\.[0-9]+","",gene),fgeneid=as.integer(gsub("ENSG","",ensid)))
predw <- mutate(predw,rsid=as.integer(gsub("rs","",rsid)))
predw <- select(predw,-gene,-ensid)
eqtldf <- left_join(eqtldf,predw,by=c("rsid","fgeneid"))
eqtldf <- mutate(eqtldf,weight=ifelse(is.na(weight),0,weight)) %>% select(-ref_allele,-eff_allele,-sgeneid,-pvalg,-seg,-betahat)

picis <- group_by(eqtldf,fgeneid) %>% summarise(nonzw=sum(weight!=0),picis=nonzw/n())
filter(picis,picis!=0) %>% ggplot()+geom_histogram(aes(x=picis),binwidth = 0.001)
eqtl_pred <- mutate(eqtl_pred,weight=ifelse(ref_allele<eff_allele,weight,-weight))
cor((eqtl_pred$weight),(eqtl_pred$thetahat))
ggplot(eqtl_pred)+stat_binhex(aes(x=weight,y=thetahat))
head(eqtl_pred)









lmv <- readRDS("~/Desktop/eQTL/lmvegas.RDS") %>% ungroup()
group_by(gwas_cis,fgeneid) %>% summarise(bestp=)








#
#
#
# serrs <- serrMatrix(osnpdat,bexpdat,betas)
# gexp <- osnpdat%*%betas
#
# geno <- SlicedData$new()
# expr <- SlicedData$new()
# cvrt <- SlicedData$new()
#
# geno$CreateFromMatrix(t(osnpdat))
# expr$CreateFromMatrix(t(bexpdat))
# #cvrt$CreateFromMatrix(t(covdat))
#
# me = Matrix_eQTL_engine(snps = geno,gene = expr,cvrt = SlicedData$new(),
#                         output_file_name = "~/Dropbox/eQTL/eQTL_test.txt",
#                         pvOutputThreshold = 1e-6,
#                         useModel = modelLINEAR,
#                         errorCovariance = numeric(),
#                         verbose = TRUE,
#                         pvalue.hist = TRUE,
#                         min.pv.by.genesnp = FALSE,
#                         noFDRsaveMemory = FALSE);
# meqtls <-me$all$eqtls
# meqtls <- mutate(meqtls,snpid=as.integer(gsub("row","",snps))-1,expid=as.integer(gsub("row","",gene))-1)
# bres <- inner_join(mdf,meqtls,by=c("snp"="snpid","gene"="expid"))


#
#
#
# ngeno <- SlicedData$new()
# nexpr <- SlicedData$new()
# ncvrt <- SlicedData$new()
#
# ngeno$CreateFromMatrix(t(ogenodat))
# nexpr$CreateFromMatrix(t(osexpdat))
#
# nme = Matrix_eQTL_engine(snps = ngeno,gene = expr,cvrt = ncvrt,
#                         output_file_name = "~/Dropbox/eQTL/neQTL_test.txt",
#                         pvOutputThreshold = 1,
#                         useModel = modelLINEAR,
#                         errorCovariance = numeric(),
#                         verbose = TRUE,
#                         pvalue.hist = TRUE,
#                         min.pv.by.genesnp = FALSE,
#                         noFDRsaveMemory = FALSE);
# nme$all$eqtls$beta-me$all$eqtls$beta



big_eqtlh5 <- "/media/nwknoblauch/Data/GTEx/all_Whole_Blood_eQTL.h5"
for(i in 1:22){
  cat(paste0(i," of 22\n"))
  trans_df <- read_h5_df(oh5s[i],"cis_eQTL")
  write_h5_df(trans_df,group="cis_eQTL",outfile = big_eqtlh5,deflate_level = 4)
}










