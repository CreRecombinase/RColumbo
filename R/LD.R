
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



gtex_gwas_intersect <- function(gtex_h5,gwas_h5s,output_h5){
  require(h5)
#  gtex_h5 <- "/media/nwknoblauch/Data/GTEx/GTEx_SNP_h5/Adipose_Subcutaneous.h5"
#  gwas_h5s <- paste0("/media/nwknoblauch/Data/1kg/1kg_",1:22,".mat")
  gtex_df <- read_df_h5(gtex_h5,"SNPinfo") %>% mutate(gtex_ind=1:n())

  gwas_df <- bind_rows(lapply(gwas_h5s,read_df_h5,subcols = c("chr","pos")))
  gtex_gwas_df <- inner_join(gtex_df,gwas_df,by=c("chrom"="chr","pos"))

  write_df_h5(gtex_gwas_df,groupname = "SNPinfo",outfile = output_h5)


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
  library(Matrix)
   # input_h5file <- "/media/nwknoblauch/Data/GTEx/1kg_SNP_H5/EUR.chr1_1kg.h5"
  # m=85
  # Ne=11490.672741
  # cutoff=1e-3
  # rowchunk <- 1
  # output_h5file <- paste0("/media/nwknoblauch/Data/GTEx/1kg_LD/EUR.chr1_",rowchunk,"_1kg.h5")
  # chunksize <- 25000


  stopifnot(file.exists(input_h5file))
  library(h5)
  inf <- h5file(input_h5file,'r')
  geno_d <- inf['/SNPdata/genotype']
  map_d <- inf['/SNPinfo/map']
  n_snps <- ncol(geno_d)
  rowchunk_tot <-ceiling(n_snps/chunksize)
  stopifnot(as.integer(chunksize/2)==chunksize/2)
  chunksize <- as.integer(chunksize/2)
  stopifnot(dim(map_d)==n_snps,rowchunk<=rowchunk_tot)
  snpAsnps <- ((rowchunk-1)*chunksize+1):min(n_snps,chunksize*rowchunk)
  stopifnot(length(snpAsnps)<=chunksize)

  snpA <- geno_d[,snpAsnps]
  Asnps <- ncol(snpA)
  mapa <- map_d[snpAsnps]
  stopifnot(length(mapa)==ncol(snpA))
  ctime <- system.time(ldA <-calcLD(hmata = snpA,hmatb = snpA,mapa = mapa,mapb = mapa,m = m,Ne = Ne,cutoff = cutoff,isDiag = T))
  ld_sp <- gen_sparsemat(ldmat = ldA,istart = snpAsnps[1],jstart = snpAsnps[1],nSNPs = n_snps)

  rm(ldA)
  gc()
  if(snpAsnps[length(snpAsnps)]<n_snps){
    snpBsnps <- (snpAsnps[length(snpAsnps)]+1):min(n_snps,chunksize*(rowchunk+1))
    stopifnot(length(snpBsnps)<=chunksize)
    snpB <- geno_d[,snpBsnps]
    mapb <- map_d[snpBsnps]
    stopifnot(length(mapb)==ncol(snpB))
    ldB <- calcLD(hmata=snpB,hmatb=snpB,mapa=mapb,mapb=mapb,m=m,Ne=Ne,cutoff=cutoff,isDiag=T)
    ldB_sp <- gen_sparsemat(ldmat = ldB,istart = snpBsnps[1],jstart = snpBsnps[1],nSNPs = n_snps)
    rm(ldB)
    gc()
    ld_sp <- ld_sp+ldB_sp
    ldAB <- calcLD(hmata=snpA,hmatb=snpB,mapa=mapa,mapb=mapb,m=m,Ne=Ne,cutoff=cutoff,isDiag=F)
    rm(snpA,snpB)
    gc()
    ldAB_sp <- gen_sparsemat(ldmat = ldAB,istart = snpAsnps[1],jstart = snpBsnps[1],nSNPs = n_snps)
    ld_sp <- ld_sp+ldAB_sp
    snpAsnps <- c(snpAsnps,snpBsnps)
  }
  h5close(inf)
  write_ccs_h5(h5filename = output_h5file,spmat = ld_sp,groupname = "R",symmetric = T)
  gc()
  return(snpAsnps)
}


concat_LD <- function(h5filenames,output_h5filename,groupname="R"){
  # output_h5filename <- "/media/nwknoblauch/Data/GTEx/1kg_LD/EUR.chr19_1kg.h5"
  tot_chunks <- length(h5filenames)
  library(h5)
  for(i in 1:tot_chunks){
    cat(i,"\n")
    cat(h5filenames[i],"\n")
    cat("i,")
    append_h5_h5(h5filenames[i],
                 input_groupname = "R",
                 input_dataname = "i",
                 output_h5 = output_h5filename,
                 output_groupname = groupname,
                 output_dataname = "i")
    cat("j,")
    append_h5_h5(h5filenames[i],
                 input_groupname = "R",
                 input_dataname = "j",
                 output_h5 = output_h5filename,
                 output_groupname = groupname,
                 output_dataname = "j")
    cat("data\n")
    append_h5_h5(h5filenames[i],
                 input_groupname = "R",
                 input_dataname = "data",
                 output_h5 = output_h5filename,
                 output_groupname = groupname,
                 output_dataname = "data")
    gc()
  }
  inf <- h5file(h5filenames[1],'r')
  ing <- inf[groupname]
  data_dim <- h5attr(ing,"Dimensions")
  layout <- h5attr(ing,"Layout")
  h5close(inf)
  of <-h5file(output_h5filename,'a')
  ofg <- of[groupname]
  h5::createAttribute(.Object = ofg,attributename = "Layout",data="COO")
  h5::createAttribute(.Object = ofg,attributename = "Dimensions",data=  data_dim)
  h5close(of)
}




