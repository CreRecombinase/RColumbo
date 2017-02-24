library(gdsfmt)
library(SNPRelate)
library(RColumbo)
library(dplyr)
library(readr)
library(tidyr)





dgn_gdsnf <- "/media/nwknoblauch/Data/DGN/Lev/DGN.gds"
dgn_h5 <- "/media/nwknoblauch/Data/DGN/Lev/eQTL/DGN_ortho.h5"
dgn_bed <- "/media/nwknoblauch/Data/DGN/Lev/DGN"
#gtex_h5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data.h5"
#gtex_bed <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed"
dgn_famf <- "/media/nwknoblauch/Data/DGN/Lev/DGN.fam"
expf <- "/media/nwknoblauch/Data/DGN/Lev/data_used_for_eqtl_study/cis_data.txt"
#gtex_famf <- "/media/nwknoblauch/Data/GTEx/Whole_Blood_v6p_bed.fam"
expdir <- "/home/nwknoblauch/Desktop/eQTL/Heritability/expdir/"
dgn_gds <- snpgdsOpen(dgn_gdsnf)
all_sample_ids <- data_frame(whole_id=read.gdsn(index.gdsn(dgn_gds,"sample.id")))
all_sample_ids <- separate(all_sample_ids,col = whole_id,into = c("first","second","third","fourth","fifth"),remove = F)


ex_ids <- read_delim(expf,delim = "\t",col_names="ID")

ex_ids <- filter(ex_ids,!is.na(ID))
sub_sample_ids <- semi_join(all_sample_ids,ex_ids,by=c("fourth"="ID"))
#gtexgds <- snpgdsOpen(gtex_gdsnf)
snpgdsClose(dgn_gds)
snpgdsGDS2BED(dgn_gdsnf,dgn_bed,sample.id = sub_sample_ids$whole_id,)

nexp <- get_rownum_h5(dgn_h5,"EXPdata","expression")
expmat <- read_fmat_h5(hap_h5file = dgn_h5,groupname = "EXPdata",dataname = "expression",offset = 0,chunksize = nexp)
sexpmat <- scale(expmat,center = T,scale = T)
expleg <- read_h5_df(dgn_h5,"EXPinfo")

expped <- read.table(dgn_famf,header=F,sep="\t")
colnames(expped) <-c("Family","Individual","Paternal","Maternal","Sex",
                     "Expression")
i <- 1

grmc <- paste0("/home/nwknoblauch/Downloads/gcta/gcta64 --bfile ",dgn_bed,
               " --maf 0.01 --make-grm --out test --thread-num 10")


setwd("/media/nwknoblauch/Data/DGN/Lev")
system(grmc)

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
  if(file.exists("ptest.hsq")){
    fdf <-read.table("ptest.hsq",sep="\t",nrows = 4,header=T)
    fdfl[[i]] <- fdf
    ndf <- read.table("ptest.hsq",sep="\t",skip = 5,header=F)
    ndfl[[i]] <- ndf
    file.remove(expf)
    file.remove("ptest.hsq")
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




