
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



flipmat <- function(hapd,doFlip){
  nhapd <- rbind(doFlip,hapd)
  rhapd <- apply(nhapd,2,function(x){
    if(x[1]==0){
      return(x[-1])
    }else{
      return(abs(1-x[-1]))
    }
  })
  return(rhapd)
}






chunk_df <- function(df,chunk.size=NULL,n.chunks=NULL){
  require(BBmisc)
  require(dplyr)
  stopifnot(xor(!is.null(chunk.size),!is.null(n.chunks)))

  ind <-1:nrow(df)
  if(!is.null(chunk.size)){
    if(chunk.size==nrow(df)){
      dfl <- list()
      dfl[[1]]<- df
      return(dfl)
    }
      ichunk <- chunk(ind,chunk.size = chunk.size)
  }else{
    if(n.chunks==1){
      dfl <- list()
      dfl[[1]]<- df
      return(dfl)
    }
    ichunk <- chunk(ind,n.chunks = n.chunks)
  }
  dfl <- list()
  for(i in 1:length(ichunk)){
    dfl[[i]] <- dplyr::slice(df,ichunk[[i]])
  }
  return(dfl)
}




chunk_eQTL_mat <- function(exph5,snph5,outh5,snpinter=NULL,expinter=NULL){
  require(h5)
  #cat("Starting to read snph5:",snph5)
  #snpleg <-read_df_h5(snph5,"SNPinfo",filtervec=snpinter)
  #cat("Starting to read exph5:",exph5)
  #expleg <- read_df_h5(exph5,"EXPinfo",filtervec=expinter)
  #expdat <- read_dmat_chunk_ind(exph5,"EXPdata","expression",expinter)
  #snpdat <- read_dmat_chunk_ind(snph5,"SNPdata","genotype",snpinter)
  #  feqtl <- really_fast_eQTL(Genotype = snpdat,snpanno = snpleg,Expression = expdat,expanno = expleg)
  cat("Reading Genotype data\n")
  genof <-h5file(snph5,mode = 'r')
  Genotype <- genof["/SNPdata/genotype"][,snpinter]
  h5close(genof)
  cat("Reading Expression\n")
  expf <- h5file(exph5,mode='r')
  Expression <- expf["/EXPdata/expression"][,expinter]
  h5close(expf)
  cat("Mapping eQTL\n")
  eqtl <- fastest_eQTL(Genotype,Expression)
  cat("Reading legends\n")
  sub_expleg <- read_df_h5(exph5,"EXPinfo",filtervec=expinter)
  sub_snpleg <- read_df_h5(snph5,"SNPinfo",filtervec=snpinter)
  cat("Writing eQTLmats\n")
  write_2dmat_h5(h5f = outh5,groupn = "eQTL",datan = "beta_mat",chunksize = as.integer(c(length(snpinter)/2,length(expinter)/2)),deflate_level = 4,data = eqtl[,,1])
  write_2dmat_h5(h5f = outh5,groupn = "eQTL",datan = "t_mat",chunksize = as.integer(c(length(snpinter)/2,length(expinter)/2)),deflate_level = 4,data = eqtl[,,2])
  cat("Writing legends\n")
  write_df_h5(df = sub_snpleg,groupname = "SNPinfo",outfile = outh5,deflate_level = 4)
  write_df_h5(df = sub_expleg,groupname = "EXPinfo",outfile = outh5,deflate_level = 4)
  cat("Done!\n")
  return(dim(eqtl))
}



block_LD <- function(input_h5file,output_h5file,m=85,Ne=11490.672741,cutoff=1e-3,rowchunk,chunksize){

  input_h5file <- "/media/nwknoblauch/Data/GTEx/1kg_SNP_H5/EUR.chr19_1kg_map.h5"
  output_h5file <- "/media/nwknoblauch/Data/GTEx/1kg_LD/EUR.chr19_1kg.h5"
  m=85
  Ne=11490.672741
  cutoff=1e-3
  rowchunk <- 1
  chunksize <- 25000
  stopifnot(file.exists(input_h5file))
  library(h5)
  inf <- h5file(input_h5file,'r')
  geno_d <- inf['/SNPdata/genotype']
  map_d <- inf['/SNPinfo/map']
  n_snps <- ncol(geno_d)
  stopifnot(dim(map_d)==n_snps)
  snpAsnps <- ((rowchunk-1)*chunksize+1):min(n_snps,chunksize*rowchunk)
  stopifnot(length(snpAsnps)<=chunksize)

  snpA <- geno_d[,snpAsnps]
  Asnps <- ncol(snpA)
  mapa <- map_d[snpAsnps]
  stopifnot(length(mapa)==ncol(snpA))
  distmat_0 <- matrix(0,nrow = Asnps,ncol = Asnps)
  ctime <- system.time(ldA <-calcLD(hmata = snpA,hmatb = snpA,mapa = mapa,mapb = mapa,m = m,Ne = Ne,cutoff = cutoff,isDiag = T))
  ldsp <- gen_sparsemat(ldmat = ldA,istart = 0,jstart = 0,nSNPs = ncol(ldA))


}





