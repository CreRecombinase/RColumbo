library(RColumbo)
library(dplyr)
library(rhdf5)
library(SQUAREM)
#IBDs <- paste0("/media/nwknoblauch/Data/IBD_GWAS_chr",1:22,".RDS")
# fxibdf <- dir("~/Desktop/eQTL/round2_squarem_step2/",full.names=T)
# readfxibdf <- function(xibdfile){
#   tdf <- data_frame(h=c(h5read(xibdfile,"h")),
#                     logw=c(h5read(xibdfile,"logw")),
#                     theta0=c(h5read(xibdfile,"theta0")))
#   return(tdf)
# }
# xibd <- bind_rows(lapply(fxibdf,readfxibdf))




hfiles <- paste0("/media/nwknoblauch/Data/rss/heighth/height2014.chr",1:22,".mat.h5")
samdf <- lapply(1:22,function(x,ibdfiles){
  cat(x,"\n")
  tdf <- data_frame(logpi=c(h5read(ibdfiles[x],name="logpisam")),
                    hsam=c(h5read(ibdfiles[x],name="hsam"))) %>% mutate(chrom=x)
  return(tdf)},ibdfiles=fxibdf)
samdf <- bind_rows(samdf)
#tgs <- h5read(fxibdf[1],name="gammasam")
samdf <- mutate()
ggplot(samdf)+geom_histogram(aes(x=logpi,..density..))+facet_wrap(~chrom,scales="free")


alpisam <- lapply(fxibdf,function(x){c(h5read(x,name="logpisam"))})
hist(unlist(alpisam))
xibd <- dir("~/Desktop/eQTL/round2_squarem_step2/",full.names = T)

xheight <- paste0("~/Desktop/eQTL/height_res_h5/height2014.chr",1:22,".beta.mat.h5")
xhdf <- bind_rows(lapply(xheight,read_h5_df,groupname=""))
xhdf <- mutate(xhdf,finid=1:n())
ppdf <- inner_join(tgs,xhdf,by="finid") %>% select(-chrom,-finid)
mutate(ppdf,Z=rbinom(n(),size=1,prob=mugam)) %>% filter(Z==1) %>% summarise(sa=sum((beta_bvsr_mean^2-se^2)/n()))
#filter(ppdf,)
ppdf <- mutate()
xhdf <- mutate(xhdf,bvsr_z=ifelse(beta_bvsr_median==0,T,F))



xibd_df <- bind_rows(lapply(xibd,function(x){
  data_frame(logw=c(h5read(x,"logw")),theta0=c(h5read(x,"theta0")))
}))


mutate(xibd_df,c=max(logw),w=exp(logw-c),nw=w/sum(w)) %>% summarise(ntheta0=10^sum(nw*theta0))

summarise(xibd_df,pi=sum(theta0*exp(logw)))

all_gwas <- bind_rows(lapply(IBDs,readRDS))
all_gwas <- mutate(all_gwas,group=as.integer(group))
xh5 <- "~/Desktop/eQTL/ibd2015_sumstat_chr19.h5"
nxh5 <- "~/Desktop/eQTL/ibd2015.chr1.mat.h5"
h5ls(nxh5)
xpos <- c(h5read(xh5,"pos_uni"))
xbeta <- c(h5read(xh5,"betahat_uni"))
xdf <- data_frame(pos=xpos,xbeta=xbeta)

x_all <- filter(all_gwas,chrom==19) %>% inner_join(xdf)
plot(x_all$beta,x_all$xbeta)
ogwasd <- read.table("/media/nwknoblauch/Data/gwas_data/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc.gz",header=T,stringsAsFactors = F)
ogwasd <- mutate(ogwasd,ntot=nchar(Direction),
                 fpos=str_count(Direction,pattern=fixed("+"))/ntot,
                 fneg=str_count(Direction,pattern=fixed("-"))/ntot,
                 fques=str_count(Direction,pattern=fixed("?"))/ntot,
                 fzero=str_count(Direction,pattern=fixed("0"))/ntot)
hist(ogwasd$fpos/ogwasd$fneg)
ggwas <- filter(ogwasd,fpos+fneg==1)
hist(ggwas$fpos)
all_gwas <- mutate(all_gwas,SNP=paste0("rs",rsid))
nall_gwas <- inner_join(all_gwas,ggwas)



emrl <- list()
nsnps <- numeric(50)
for(i in 1:50){
  cat(paste0("i:",i,"\n"))
  sub_snp <- group_by(nall_gwas,chrom,group) %>% sample_n(1) %>% ungroup()
  nsnps[i] <- nrow(sub_snp)
  emrl[[i]] <- list()
  p <- nrow(sub_snp)
  for(j in 1:10){
  # cat(paste0("j:",j,"\n"))
    npi <- c(exp(runif(1,log(1/p),0)),exp(runif(1)))
    emrl[[i]][[j]] <- squarem(par=npi,bh=sub_snp$beta,si=sub_snp$serr,fixptfn=sslab_em)
  }
}

pimat <- matrix(0,50,10)
taumat <- matrix(0,50,10)



lpi <- exp(runif(10000,log(1/p),0))
opi <- runif(10000,1/p,1)
nopi <- runif(10000,1/p,1)
qqplot(opi,lpi)
for(i in 1:50){
  for(j in 1:10){
    pimat[i,j] <- emrl[[i]][[j]][["par"]][1]
    taumat[i,j] <- emrl[[i]][[j]][["par"]][2]
  }
}

plot(c(pimat),c(taumat))
hist(c(taumat))
hist(c(pimat))
hist(nsnps)
