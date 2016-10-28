#Code to simulate GWAS data

betasim <- function(annodf){
  require(dplyr)
  resdf <- mutate(annodf,Z=rbinom(n(),size = 1,prob=pi)) %>% mutate(Beta=ifelse(Z==1,rnorm(n=n(),mean=0,sd=sa),0))
  return(resdf)
}

betasim_h5 <- function(h5file,pi,sa){
  legdf <- read_h5_df(h5file,"Legend") %>% mutate(betahat=beta)
  resdf <- mutate(legdf,Z=rbinom(n(),size = 1,prob=pi)) %>% mutate(Beta=ifelse(Z==1,rnorm(n=n(),mean=0,sd=sa),0))
  return(resdf)
}
#

yhat <-function(tdf){
  #chunksize <- get_rownum_h5(tdf$h5file[1],"eQTL","genotype")
  #tgeno <- read_fmat_h5(tdf$h5file[1],"eQTL","genotype",0,chunksize)
  tgeno <- read_fmat_chunk_ind(tdf$h5file[1],"eQTL","genotype",tdf$ind)
  tgeno <-scale(tgeno,center = T,scale = F)
  ty <- tgeno %*%  tdf$Beta
  retdf <- data_frame(ind=1:length(ty),yh=c(ty),chrom=tdf$chrom[1])
  return(retdf)
}


txbeta <- function(tdf){
  # cat(tdf$chrom[1])
  tgeno <- read_geno_h5_pos(tdf$h5file[1],tdf)
  tgeno <- scale(tgeno,center=T,scale=F)
  ntdf <- filter(tdf,pos %in% as.integer(colnames(tgeno)))
  tgeno <- tgeno[,as.character(ntdf$pos)]
  ty <- tgeno %*% ntdf$Beta
  retdf <- data_frame(ind=1:length(ty),yh=c(ty),chrom=ntdf$chrom[1])
  return(retdf)
}


xbetasim <- function(betadf,h5file){
 #h5filedf <-data_frame(chrom=1:length(h5filevec),h5file=h5filevec)
  # betadf <- mutate(betadf,h5file=h5file)
  tydf <- group_by(betadf,chrom) %>% do(yhat(.))
  fy <- group_by(tydf,ind) %>% summarise(xb=sum(yh))
  return(fy)
}

doGWAS <- function(h5files,pi,sa,chunksize=100000){
  stopifnot(length(h5files)==22)
  betadf <- bind_rows(lapply(h5files,betasim_h5,pi=pi,sa=sa))
  betadf <- group_by(betadf,chrom) %>% mutate(h5file=h5files[chrom]) %>% ungroup()
  yhat <- group_by(betadf,chrom) %>% do(yhat(.))
  ny <- group_by(yhat,ind) %>% summarise(xb=sum(yh))

  fb <- list()
  for(i in 1:length(h5files)){
    cat(i,"\n")
    fb[[i]] <- ffast_GWAS(h5files[i],phenotype = ny$xb,chunksize = chunksize) %>% mutate(chrom=i)
  }
  fbdf <- bind_rows(fb)
  retdf <- select(betadf,chrom,pos,Z,Beta) %>% mutate(betahat=fbdf$Betahat)
  # gbeta <- filter(retdf,Z!=0)
  return(retdf)
}
# group_by(retdf,chrom) %>% summarise(nsnp=n()) %>% arrange(desc(nsnp))






