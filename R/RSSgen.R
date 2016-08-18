gen_RSSmat <- function(chunkRDSf,gwash5,outh5,Nind=34652){
  require(Matrix)
  require(rhdf5)
  require(dplyr)

  bmat <-readRDS(chunkRDSf)
  bmat<- as(bmat,"dgTMatrix")
  ind <- unique(c(bmat@i,bmat@j))
  bmat <- bmat[ind+1,ind+1]
  bmat <- as.matrix(bmat)
  legdf <- read_h5_df(haph5,"Legend")
  gwasd <- read_h5_df(gwash5,"GWAS")
  gwasd <- mutate(gwasd,rsid=paste0("rs",rsid))

  rslistvec <- subset_h5_ref(haph5,mapfile = mapfile,eqtlfile)
  readind <-which((legdf$rsid %in% rslistvec)&(!duplicated(legdf$rsid)))
  slegdf <- slice(legdf,readind)
  slegdf <- mutate(slegdf,i=0:(length(rsid)-1))
  slegdf <- filter(slegdf,i %in% ind)
  gwasd <- semi_join(gwasd,slegdf)
  bmat <- bmat+t(bmat)
  diag(bmat) <- 1

  h5createFile(outh5)
  h5createDataset(outh5,"R_uni",dims=dim(bmat),storage.mode="double",chunk=c(1000,1000),level=2)
  h5createDataset(outh5,"betahat_uni",dims=nrow(gwasd),storage.mode="double",chunk=c(1000),level=2)
  h5createDataset(outh5,"se_uni",dims=nrow(gwasd),storage.mode="double",chunk=c(1000),level=2)
  h5createDataset(outh5,"Nsnp",dims=nrow(gwasd),storage.mode="double",chunk=c(1000),level=2)
  h5createDataset(outh5,"rsid",dims=nrow(gwasd),storage.mode="character",size=max(nchar(gwasd$rsid))+1,chunk=c(1000),level=2)

  gwasd <- mutate(gwasd,Nsnp=Nind)
  h5write(gwasd$Nsnp,file=outh5,name="Nsnp")
  h5write(gwasd$beta,file=outh5,name="betahat_uni")
  h5write(gwasd$serr,file=outh5,name="se_uni")
  h5write(gwasd$rsid,file=outh5,name="rsid")
  h5write(obj=bmat,file=outh5,name="R_uni")

}
