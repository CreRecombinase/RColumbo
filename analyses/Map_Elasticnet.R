library(glmnet)
library(broom)
library(dplyr)
library(RColumbo)
library(BBmisc)
nchunks <- 100
mychunk <- as.integer(commandArgs(trailingOnly = F))
workdir <- "/home/nwknoblauch/Desktop/eQTL/Snake/"
setwd(workdir)
rawh5 <- "Whole_Blood_eQTL_v6p_raw_data.h5"
cistf <- "GTEx_cistrans.h5"
stopifnot(file.exists(rawh5),file.exists(cistf))
cisleg <- read_h5_df(cistf,"cis")
expleg <- distinct(cisleg,exp_ind,.keep_all = T)
expdat <- read_fmat_chunk_ind(h5file = rawh5,groupname = "EXPdata","orthoexpression",expleg$exp_ind)
colnames(expdat) <- as.character(expleg$fgeneid)
ufgid <- expleg$fgeneid
ensid <- paste0("ENSG",sprintf("%011d",ufgid))
fgidchunk <- chunk(ufgid,n.chunks = nchunks)
h5chunks <- paste0("/media/nwknoblauch/Data/GTEx/cis_elasticnet/chunk_",1:length(fgidchunk),".h5")
mh5 <- h5chunks[mychunk]
mufgid <- fgidchunk[[mychunk]]
ensid <- paste0("ENSG",sprintf("%011d",mufgid))

overall_tl <- list()
lambda_min_tl <- list()
weights_tl <- list()
bcisleg <- filter(cisleg,fgeneid %in% mufgid)
rm(cisleg)
gc()
bsnpleg <- distinct(bcisleg,snp_ind,.keep_all = T) %>% mutate(snp_name=paste0("chr",chrom,"_",pos))
bsnpdat <- read_fmat_chunk_ind(rawh5,"SNPdata","genotype",bsnpleg$snp_ind)
colnames(bsnpdat) <- paste0("chr",bsnpleg$chrom,"_",bsnpleg$pos)
bcisleg <- mutate(bcisleg,snp_name=paste0("chr",chrom,"_",bcisleg$pos))
namel <- split(bcisleg$snp_name,bcisleg$fgeneid)
for(i in i:length(namel)){
  cat(ensid[i],"_",i,"\n")
  texp <- expdat[,as.character(mufgid[i])]
  snp_data <- bsnpdat[,namel[[i]],drop=F]
  if(ncol(snp_data)>1){
  pred_fit <- cv.glmnet(x = snp_data,y = texp,alpha=0.95,parallel = F)
  overall_tl[[as.character(mufgid[i])]] <- tidy(pred_fit) %>% mutate(fgeneid=mufgid[i])
  lambda_min_tl[[as.character(mufgid[i])]] <- glance(pred_fit) %>% mutate(fgeneid=mufgid[i])
  weights_tl[[as.character(mufgid[i])]] <- tidy(pred_fit$glmnet.fit) %>% mutate(fgeneid=mufgid[i]) %>% filter(estimate!=0)
  }else{
    overall_tl[[as.character(mufgid[i])]] <- NULL
    lambda_min_tl[[as.character(mufgid[i])]] <- NULL
    weights_tl[[as.character(mufgid[i])]] <- NULL
  }
}
overall <- bind_rows(overall_tl)
lambda_min <- bind_rows(lambda_min_tl)
weights <- select(bsnpleg,chrom,pos,snp_name) %>% right_join(bind_rows(weights_tl),by=c("snp_name"="term"))
weights <- mutate(weights,chrom=ifelse(is.na(chrom),0,chrom),pos=ifelse(is.na(pos),0,pos)) %>% select(-snp_name)

write_h5_df(weights,group = "weights",outfile = mh5,deflate_level = 4)
write_h5_df(overall,group = "overall",outfile = mh5,deflate_level = 4)
write_h5_df(lambda_min,group = "lambda_min",outfile = mh5,deflate_level = 4)



