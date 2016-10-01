library(RColumbo)
library(dplyr)
library(rhdf5)
library(SQUAREM)
library(readr)
#IBDs <- paste0("/media/nwknoblauch/Data/IBD_GWAS_chr",1:22,".RDS")
# fxibdf <- dir("~/Desktop/eQTL/round2_squarem_step2/",full.names=T)
# readfxibdf <- function(xibdfile){
#   tdf <- data_frame(h=c(h5read(xibdfile,"h")),
#                     logw=c(h5read(xibdfile,"logw")),
#                     theta0=c(h5read(xibdfile,"theta0")))
#   return(tdf)
# }
# xibd <- bind_rows(lapply(fxibdf,readfxibdf))

genofile <- "/media/nwknoblauch/Data/DGN/Lev/testGenRed.II.autosomal.ped.gen.gz"
mapfile <- "/media/nwknoblauch/Data/DGN/Lev/testGenRed.II.autosomal.map"
dbsnpfile <- "/home/nwknoblauch/Desktop/eQTL/Snake/dbsnp.h5"
Nind <- 2226
Nsnp <- 720591
dgn_h5 <- "/home/nwknoblauch/Desktop/eQTL/DGN.h5"
haplotype_h5s <-paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
thaph5 <- haplotype_h5s[1]
h5ls(thaph5)
#write_genofile_h5(genofile = genofile,dbsnpfile = dbsnpfile,Nind = 2226,Nsnps = Nsnp,chunksize = 10000,h5file = dgn_h5,deflate_level = 2)

height_raw_files <- paste0("/media/nwknoblauch/Data/rss/heighth/height2014.chr",1:22,".mat.h5")
height_res_files <- paste0("/media/nwknoblauch/Data/rss/height/rssbvsr_height2014_chr",1:22,"_result.mat")

#betasam <- colMeans(h5read(height_res_files[1],name="betasam"))

posdf <- bind_rows(lapply(1:22,function(x,ibdfiles){
  cat(x,"\n")
  tdf <- data_frame(chrom=c(h5read(ibdfiles[x],name="chr")),
                    pos=c(h5read(ibdfiles[x],name="pos19")))
  return(tdf)
},ibdfiles=height_raw_files))


samdf <- bind_rows(lapply(1:22,function(x,ibdfiles){
  cat(x,"\n")
  tdf <- data_frame(pi=mean(exp(c(h5read(ibdfiles[x],name="logpisam")))),
                    h=mean(c(h5read(ibdfiles[x],name="hsam")))) %>% mutate(chrom=x)
  return(tdf)},ibdfiles=height_res_files))

sadf <- bind_rows(lapply(1:22,function(x,ibdfiles){
  cat(x,"\n")
  tdf <- data_frame(betahat=c(h5read(ibdfiles[x],name="betahat")),
                      se=c(h5read(ibdfiles[x],name="se")),
                    N=c(h5read(ibdfiles[x],name="Nsnp"))) %>% mutate(chrom=x,s=sqrt(se^2+betahat^2/N))
  return(tdf)},ibdfiles=height_raw_files))
sedf <- inner_join(samdf,sadf)
hedf <- group_by(sedf,chrom) %>% summarise(h=h[1],pi=pi[1],sa=sqrt(h/(pi*sum(1/(N*s^2)))))
fiparam <-summarise(sedf,sa=mean(sqrt(h/(pi*sum(1/(N*s^2))))),pi=mean(pi))


h5files <- paste0("/media/nwknoblauch/Data/GTEx/chr",1:22,"_1kg_Whole_Blood_HEIGHT.h5")
tbdf <- bind_rows(lapply(h5files,betasim_h5,pi=pi,sa=sa))
tbdf <- group_by(tbdf,chrom) %>% mutate(h5file=h5files[chrom]) %>% ungroup()

tbdf <- betasim_h5(h5file = th5file,pi = pi,sa = sa)

annodf <- inner_join(hedf,posdf)
annodf <- mutate(annodf,pi=mean(pi),sa=mean(sa))
betadf <- betasim(annodf)

geno_df <- read_h5_df(dgn_h5,"Legend") %>% mutate(index=1:n()) %>%select(-rsid)
nbetadf <- inner_join(geno_df,betadf) %>% distinct(chrom,pos,.keep_all=T)
ny <- xbetasim(nbetadf,dgn_h5)
indvec <- nbetadf$index

no_err <- fast_GWAS(indvec,dgn_h5,ny$xb,chunksize=100000)



