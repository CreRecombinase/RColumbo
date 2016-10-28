library(RColumbo)
library(dplyr)
library(rhdf5)
library(SQUAREM)
library(readr)
library(ggplot2)
#IBDs <- paste0("/media/nwknoblauch/Data/IBD_GWAS_chr",1:22,".RDS")
# fxibdf <- dir("~/Desktop/eQTL/round2_squarem_step2/",full.names=T)
# readfxibdf <- function(xibdfile){
#   tdf <- data_frame(h=c(h5read(xibdfile,"h")),
#                     logw=c(h5read(xibdfile,"logw")),
#                     theta0=c(h5read(xibdfile,"theta0")))
#   return(tdf)
# }
# xibd <- bind_rows(lapply(fxibdf,readfxibdf))
oh5files <- paste0("/media/nwknoblauch/Data/DGN/LDmats/chr",1:22,"_1kg_DGN_Height_LD.h5")
# genofile <- "/media/nwknoblauch/Data/DGN/Lev/testGenRed.II.autosomal.ped.gen.gz"
# mapfile <- "/media/nwknoblauch/Data/DGN/Lev/testGenRed.II.autosomal.map"
# dbsnpfile <- "/home/nwknoblauch/Desktop/eQTL/Snake/dbsnp.h5"
# Nind <- 2226
# Nsnp <- 720591

dgn_h5 <- "/home/nwknoblauch/Desktop/eQTL/DGN.h5"
# haplotype_h5s <-paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
# thaph5 <- haplotype_h5s[1]
# h5ls(thaph5)
#write_genofile_h5(genofile = genofile,dbsnpfile = dbsnpfile,Nind = 2226,Nsnps = Nsnp,chunksize = 10000,h5file = dgn_h5,deflate_level = 2)

height_raw_files <- paste0("/media/nwknoblauch/Data/rss/heighth/height2014.chr",1:22,".mat.h5")
height_res_files <- paste0("/media/nwknoblauch/Data/rss/height/rssbvsr_height2014_chr",1:22,"_result.mat")

#betasam <- colMeans(h5read(height_res_files[1],name="betasam"))

sedf <- inner_join(samdf,sadf)
hedf <- group_by(sedf,chrom) %>% summarise(h=h[1],pi=pi[1],sa=sqrt(h/(pi*sum(1/(N*s^2)))))
fiparam <-summarise(sedf,sa=mean(sqrt(h/(pi*sum(1/(N*s^2))))),pi=mean(pi))
sa <- fiparam$sa
pi <- fiparam$pi

plot(hedf$pi,hedf$sa)
abline(h=sa)
abline(v=pi)

pi <- pi*2
h5files <- paste0("/media/nwknoblauch/Data/DGN/EUR.chr",1:22,"_1kg_DGN_Height.h5")
gwasdf <- bind_rows(lapply(h5files,read_h5_df,groupname="Legend"))

colordf <- mapply(function(x,chro){
  retdf <- read_h5_df(x,groupname = "Color") %>% select(-ind) %>% mutate(chrom=chro)
  return(retdf)
},oh5files,1:22,SIMPLIFY = F)
colordf <- bind_rows(colordf)
colordf <- rename(colordf,ind=index) %>% mutate(ind=ind+1)
gwasdf <- inner_join(colordf,gwasdf)

gpi <- pi
gsa <- sa
p <-c(gpi,gsa)
bh <- gwasdf$beta
si <- gwasdf$serr
#fbagwas <- mutate(fbagwas,mu=pMu(Betahat,serr = serr,pi = gpi,tau = gsa))
npar <- squarem(par = p,fixptfn = sslab_em,bh=bh,si=si)$par






pi <-exp(runif(50,log(1/p),0))


for(iter in 1:50){
#  tbdf <- bind_rows(lapply(h5files,betasim_h5,pi=pi,sa=sa))
  icdf <- filter(colordf,color==0) %>% mutate(pi=pi,sa=sa)
  icdf <- betasim(icdf) %>% mutate(ind=ind+1)
  tbdf <- icdf
  # filter(tbdf,Z==1) %>% summarise(sd(Beta),mean(Beta),sa=sa)
  tbdf <- group_by(tbdf,chrom) %>% mutate(h5file=h5files[chrom]) %>% ungroup()
  byh <- group_by(tbdf,chrom) %>% do(yhat(.))
  fyh <- ungroup(byh) %>% group_by(ind) %>% summarise(nyh=sum(yh))
  fyh <- mutate(fyh,nyh=nyh+rnorm(n(),mean=0,sd=sa/2))


  chromdf <- data_frame(chrom=1:22)
  fdf <- mutate(fyh,indc=NA) %>% full_join(mutate(chromdf,indc=NA)) %>% select(-indc) %>% group_by(chrom) %>% mutate(h5file=h5files[chrom]) %>% ungroup()
  agwas <- group_by(fdf,chrom) %>% do(ffast_GWAS(h5file=.$h5file[1],phenotype=.$nyh,chunksize=10000)) %>% ungroup()
  iagwas <- group_by(fdf,chrom)
   agwas <- group_by(agwas,chrom) %>% mutate(ind=(1:n())-1) %>% ungroup()
   bagwas <- inner_join(agwas,colordf)


filter(nagwas,Z==1) %>% summarise(sd(Betahat))

   fbagwas <- filter(bagwas,color==0)
   gpi <- pi
   gsa <- sa
   p <-c(gpi,gsa)
   bh <- fbagwas$Betahat
   si <- fbagwas$serr
   fbagwas <- mutate(fbagwas,mu=pMu(Betahat,serr = serr,pi = gpi,tau = gsa))
   npar <- squarem(par = p,fixptfn = sslab_em,bh=agwas$Betahat,si=agwas$serr)$par
nagwas <- inner_join(tbdf,agwas,by=)
  squarem(par = p,fixptfn = sslab_em,bh=fbagwas$Betahat,si=fbagwas$serr)
  ggplot(nagwas)+geom_point(aes(x=abs(Betahat)/serr,y=mu,col=as.factor(Z)))+facet_wrap(~Z)
filter(nagwas,Z==1)  %>%  ggplot(aes(x=Betahat,y=Beta))+geom_point()+geom_abline()



pivec <- numeric(100)
tauvec <- numeric(100)
pivec[1] <- pi
tauvec[1] <- sa
for(i in 2:length(pivec)){
  ret <- sslab_em(c(pivec[i-1],tauvec[i-1]),bh,si)
  pivec[i] <- ret[1]
  tauvec[i] <- ret[2]
}
plot(pivec)
plot(tauvec)

fpiter(par=p,fixptfn=sslab_em,bh=bh,si=si)
pmu <-pMu(betahat = bh,serr = si,pi = npar[1],tau = npar[2])
fbagwas <- mutate(fbagwas,mu=pmu)
fbagwas <- mutate(fbagwas,mmu=(pi*dnorm(Betahat,0,sqrt(sa^2+si^2)))/((pi*dnorm(Betahat,0,sqrt(sa^2+si^2)))+(1-pi)*dnorm(Betahat,0,si)))
ggplot(fbagwas)+geom_point(aes(x=log(percent_rank(abs(Betahat)/serr)),y=mu))
ggplot(fbagwas)+geom_histogram(aes(x=abs(Betahat)/serr))
ggplot(fbagwas)+geom_histogram(aes(x=log(mu)))

sslab_em(p,fbagwas$Betahat,si=fbagwas$serr)

sslab.lik <- function(p,bh,si){
  temp_pi <- p[1]
  temp_tau <- p[2]
  uzin <- temp_pi*dnorm(bh,mean=0,sd=sqrt(temp_tau^2+si^2))
  uizd <- uzin+(1-temp_pi)*dnorm(bh,mean=0,sd=si,log=F)
  uiz <- uzin/uizd
  fiprobs <- uiz*dnorm(bh,mean=0,sd=sqrt(temp_tau^2+si^2))
  siprobs <-(1-uiz)*dnorm(bh,mean=0,sd=si)
  iprobs <- uiz*log(dnorm(bh,mean=0,sd=sqrt(temp_tau^2+si^2))*temp_pi)+(1-uiz)*log(dnorm(bh,mean=0,sd=si)*(1-temp_pi))
  return(-sum(iprobs))
}

optim(par=p,fn=sslab.lik,bh=bh,si=si,lower = c(1e-5,1e-7),upper = c(1-1e-5,100),method="L-BFGS-B",control=list(trace=T))

pis <- seq(1e-5,1-1e-5,length.out = 20)
taus <- sort(exp(runif(20)))

lmat <- matrix(0,20,20)
for(i in 1:length(pis)){
  for(j in 1:length(taus)){
    lmat[i,j] <- sslab.lik(c(pis[i],taus[j]),bh,sqrt(si))
  }
}

betahat <- runif(100)
si <- runif(100)
mu <-runif(100)
tau <- runif(1)

tfun <- function(tau,betahat,si,mu){
  sum(-(mu*tau)/(si^2+tau^2)+(betahat^2*mu)/((si^2+tau^2)^2))
}
tfun(tau,betahat,si,mu)
plot(pis,colMeans(lmat))
plot(taus,colMeans(lmat))




