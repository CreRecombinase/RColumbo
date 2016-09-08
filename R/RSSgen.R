

gen_RSSmat <- function(df,haph5,outh5){
  require(Matrix)
  require(rhdf5)
  require(dplyr)


  legdf <- read_h5_df(haph5,"Legend")

  rsidvec <- df$rsid
  fmat <- gen_LD_rsid(rsidvec,haph5)


  stopifnot(nrow(fmat)==nrow(df))
  h5createFile(outh5)
  h5createDataset(outh5,"R_uni",dims=dim(fmat),storage.mode="double",chunk=c(1000,1000),level=2)
  h5createDataset(outh5,"betahat_uni",dims=nrow(df),storage.mode="double",chunk=c(1000),level=2)
  h5createDataset(outh5,"se_uni",dims=nrow(df),storage.mode="double",chunk=c(1000),level=2)
  h5createDataset(outh5,"Nsnp",dims=nrow(df),storage.mode="double",chunk=c(1000),level=2)
  h5createDataset(outh5,"rsid",dims=nrow(df),storage.mode="character",size=max(nchar(df$rsid))+1,chunk=c(1000),level=2)
  h5write(df$Nind,file=outh5,name="Nsnp")
  h5write(df$betahat,file=outh5,name="betahat_uni")
  h5write(df$serr,file=outh5,name="se_uni")
  h5write(df$rsid,file=outh5,name="rsid")
  h5write(obj=fmat,file=outh5,name="R_uni")
}



Gen_mixnorm <- function(n,pi,tau){
  Z <- sample(0:1,size = n,replace=T,prob = c(1-pi,pi))
  si <- runif(n)
  B <- ifelse(Z,rnorm(n=n,mean = 0,sd = tau),0)
  bhat <- rnorm(n = n,mean = B,sd = si)
  return(data.frame(betahat=bhat,si=si,Z=Z))
}



