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
grm_pref <-args[3]
cat_covf <- args[4]
cont_covf <- args[5]
bed_pref <- args[6]
estfile <- args[7]
sumfile <- args[8]


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
expmat <- read_fmat_h5(hap_h5file = h5file,groupname = "EXPdata",dataname = "expression",offset = 0,chunksize = nexp)

sexpmat <- scale(expmat,center = T,scale = T)
expleg <- read_h5_df(h5file,"EXPinfo")

expped <- read.table(bed_famf,header=F,sep="\t")
colnames(expped) <-c("Family","Individual","Paternal","Maternal","Sex",
                     "Expression")
fdfl <- list()
ndfl <- list()
expf <- tempfile()
hsqf <- tempfile()
cat("Starting!\n")
for(i in 1:nexp){
  if(i%%1000==0){
    cat(i," of ",nexp,"\n")
  }

  expped <- mutate(expped,Expression=sexpmat[,i]) %>% select(Family,Individual,Expression)
  write.table(expped,file = expf,sep="\t",col.names = F,row.names = F,
              quote=F)

  gcta_command <- paste0("/home/nwknoblauch/Downloads/gcta/gcta64 --grm ",grm_pref,
                         " --pheno ",expf,
                         " --reml ",
                         " --covar ",cat_covf,
                         " --qcovar ",cont_covf,
                         " --out ",hsqf,"  --thread-num 10 --reml-maxit 200")
  system(gcta_command,ignore.stdout = T,ignore.stderr = T)
  hsqfile <- paste0(hsqf,".hsq")
  if(file.exists(hsqfile)){
    fdf <-read.table(hsqfile,sep="\t",nrows = 4,header=T)
    fdfl[[i]] <- fdf
    ndf <- read.table(hsqfile,sep="\t",skip = 5,header=F)
    ndfl[[i]] <- ndf
    file.remove(hsqfile)
    file.remove(expf)
  }else{
    cat("ERROR!!!\n")
    fdfl[[i]] <- NULL
    ndfl[[i]] <- NULL
    file.remove(expf)
  }
}


nfdfl <- mapply(function(df,gid){
  if(!is.null(df)){
    return(mutate(df,fgeneid=gid))
  }else{
    return(NULL)
  }
},fdfl,expleg$fgeneid,SIMPLIFY = F)

nndfl <- mapply(function(df,gid){
  if(!is.null(df)){
    return(mutate(df,fgeneid=gid))
  }else{
    return(NULL)
  }
},ndfl,expleg$fgeneid,SIMPLIFY = F)

vpdf <- bind_rows(nfdfl)
vpdf <- filter(vpdf,Source=="V(G)/Vp") %>% rename(h=Variance) %>% select(-Source)
sdf <- bind_rows(nndfl)
saveRDS(sdf,sumfile)
saveRDS(vpdf,estfile)





