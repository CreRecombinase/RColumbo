
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#'
#' #' Subsets reference (haplotype) panel data by either rsid or position
#' #'
#' #' @param rsids  character vector containing rsids that will be used for generating the LD matrix
#' #' @param chrom character indicating which chromosome will be used
#' #' @param positions integer vector indicating the genomic coordinates that will be used for generating the LD matrix
#' #' @param legendfile file containing legend for haplotype data
#' #' @param hapfile file containing actual haplotype data
#' #' @param mapfile file containing physical map info
subset_ref_panel <- function(rsids=character(0),mapfile,haph5){
  require(dplyr)
  require(rhdf5)
  stopifnot(file.exists(haph5),file.exists(mapfile),
            length(rsids)>0)
  mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
  legdf <- read_h5_df(haph5,"Legend")
  nSNPs <- nrow(legdf)
  rslistvec <-rsids
  rslistvec <-intersect(rslistvec,legdf$rsid)
  rslistvec <- intersect(rslistvec,mapdf$rsid)
  stopifnot(length(rslistvec)>0)
  readind = which((legdf$rsid %in%rslistvec)&!duplicated(legdf$rsid))

  hapd <- read_haplotype_ind_h5(hap_h5file = haph5,indexes = readind)
  legdf <- slice(legdf,readind)
  mapdf <- filter(mapdf,rsid %in% rslistvec)
  colnames(hapd) <-legdf$rsid
  nhapd <- rbind(as.integer(legdf$allele0<legdf$allele1),hapd)
  rhapd <- apply(nhapd,2,function(x){
    if(x[1]==1){
      return(x[-1])
    }else{
      return(abs(1-x[-1]))
    }
  })
  retlist <-list(H=rhapd,cummap=mapdf$cummap)
  return(retlist)
}

#'
#' #' Compute the intersection of genetic map,hdf5 file, and haplotype data for a particular chromosome
#' #' @param legendfile IMPUTE format legend file
#' #' @param hapfile IMPUTE format genotype file
#' #' @param mapfile file containing physical map info
#' #' @param h5file HDF5 formatted data with merged eQTL and GWAS data
subset_h5_ref <- function(legendfile,mapfile,eqtlfile){
  require(readr)
  snpvec <- h5vec(eqtlfile,"SNP","rsid")
  snpdf <- data_frame(rsid=paste0("rs",snpvec))
  lsnpvec <- h5vec(haph5,"Legend","rsid")
  legenddf <- data_frame(rsid=lsnpvec)
  #  hapmat <- data.matrix(read_delim(hapfile,delim = " ",col_names = F))
  mapdat <- read_delim(mapfile,delim=" ",col_names=c("rsid","pos","dist"))
  sub_snp <-  semi_join(snpdf,legenddf) %>%semi_join(mapdat)
  rsl <- unique(sub_snp$rsid)
  return(rsl)
}
#'
#' #' Compute the intersection of genetic map,hdf5 file, and haplotype data for a particular chromosome
#' #' @param legendfile IMPUTE format legend file
#' #' @param hapfile IMPUTE format genotype file
#' #' @param mapfile file containing physical map info
#' #' @param h5file HDF5 formatted data with merged eQTL and GWAS data
#' #' @param gfile gencode file
#' #' @param gene character specifying which gene we want the SNPs for
subset_h5_ref_gene <- function(legendfile,mapfile,h5file,gfile,gene){
  require(readr)

  # eqtlc <-group_by(eqtldf,Gene)%>% summarise(neqtl=n(),prode=sum(log(pvale),na.rm = T),prodg=sum(log(pvalg),na.rm=T))
  # eqtlc <- mutate(eqtlc,hash=as.integer(neqtl>=100))%>%arrange(desc(hash),prode)
  eqtldf <- filter(eqtldf,Gene==gene)

  legenddf <- read_delim(legendfile,delim = " ",col_names = T)
  #  hapmat <- data.matrix(read_delim(hapfile,delim = " ",col_names = F))
  mapdat <- read_delim(mapfile,delim=" ",col_names=c("rsid","pos","dist"))
  sub_snp <- semi_join(snpdf,eqtldf) %>%mutate(rsid=paste0("rs",rsid)) %>%
    semi_join(legenddf,by=c("rsid"="ID")) %>%semi_join(mapdat,by=c("rsid"="rsid"))
  rsl <- unique(sub_snp$rsid)
  return(rsl)
}
#'
#'
#' #' Compute the intersection of genetic map,hdf5 file, and haplotype data for a particular chromosome
#' #' @param rslist list of rsids to pull
#' #' @param h5file HDF5 formatted data with merged eQTL and GWAS data
#' #' @param gfile gencode file
#' #' @param gene character specifying which gene we want the SNPs for
#' subset_h5_gene <- function(rslist,h5file,gfile,gene)
#'   h5l <- h5dfl(h5file)
#'   snpdf <- h5l[["SNP"]]
#'   eqtldf <- h5l[["eQTL"]]
#'   rm(h5l)
#'   eqtldf <- mutate(eqtldf,rsid=paste0("rs",rsid)) %>% filter(rsid %in% rslist)
#'   eqtldf <- gencode_eqtl(eqtldf,gfile)
#'
#'
#' #' Turns matrix of diagonals in to sparse banded matrix
#' #'
#' #' @param bmat matrix of diagonals
#' #' @param bwd band width
#' from_band <- function(bmat,bwd){
#'   tmat <- bandSparse(n = ncol(bmat),m = ncol(bmat),k=-c(0:bwd),diagonals = t(bmat),symmetric = T)
#'   return(tmat)
#' }
#'
#'
#'
#' #' Compute sub covariance matrix and return it as a vector
#' #'
#' #' @param sH scaled version of original matrix (mean subtracted from every column)
#' #' @param i row chunk
#' #' @param j column chunk
#' #' @param chunksize size of chunk (we're lucky that the covariance matrix is square)
#' chunk_covar <- function(sHr,sHc,isDiag=F){
#'   ret <- crossprod(sHr,sHc)/(nrow(sHr)-1)
#'   if(isDiag){
#'     return(ret[upper.tri(ret)])
#'   }else{
#'     return(c(ret))
#'   }
#' }
#'
#'
#' #' Compute which rows are going to be subset
#' #'
#' #' @param i integer designation of row chunk
#' #' @param chunksize integer size of chunks
#' #' @ncols total number of columns in the  (final) covariance matrix
#' #'
subrows <- function(i,chunksize,ncols){
  rows <-((i-1)*chunksize+1):min(ncols,(i)*chunksize)
  return(rows)
}
#'
#' #' Compute which columns are going to be subset
#' #'
#' #' @param j integer designation of column chunk
#' #' @param chunksize integer size of chunks
#' #' @ncols total number of columns in the  (final) covariance matrix
#' #'
subcols <-function(j,chunksize,ncols){
  cols <-((j-1)*chunksize+1):min(ncols,(j)*chunksize)
  return(cols)
}

#'
#'
#' #' Compute chunks of the genetic map distance matrix and return them as a vector
#' #'
#' #' @param distvec orginal cumulative distance map
#' #' @param i row chunk
#' #' @param j column chunk
#' #' @param chunksize size of chunk (we're lucky that the covariance matrix is square)
#' chunk_dist <- function(distvecr,distvecc,isDiag=F){
#'   ret <- outer(distvecr,distvecc,"-")
#'   if(isDiag){
#'     return(c(ret[upper.tri(ret)]))
#'   }else{
#'     return(c(ret))
#'   }
#' }
#'
#' fast_nS <-function(H,cummap,m,Ne,cutoff,nchunks){
#'   plan(multiprocess)

#'   stopifnot(length(ivec)==length(jvec),sum(table(ivec))==sum(table(jvec)))
#'   nmsum <- sum(1/1:(2*m-1))
#'   theta <- (1/nmsum)/(2*m+1/nmsum)


#'     sHr <- H[,subrows(ivec[i],chunksize,nSNPs),drop=F]
#'     sHc <- H[,,drop=F]
#'     distvecr <-cummap[subrows(ivec[i],chunksize,nSNPs)]
#'     distvecc <- cummap[subcols(jvec[i],chunksize,nSNPs)]
#'     tl[[i]] <-future({
#'       ivec
#'       jvec
#'       Svec <- chunk_covar(sHr,
#'                   sHc,
#'                   ivec[i]==jvec[i])
#'       distvec <-chunk_dist(distvecr,distvecc,isDiag = ivec[i]==jvec[i])
#'       shrinkage <- 4*Ne*(-distvec)/100
#'       shrinkage <- exp(-shrinkage/(2*m))
#'       shrinkage[shrinkage<cutoff] <- 0
#'       nS <- shrinkage*Svec
#'       rows <- subrows(ivec[i],chunksize = chunksize,ncols = nSNPs)
#'       cols <- subcols(jvec[i],chunksize = chunksize,ncols = nSNPs)
#'       tmat <- matrix(0,length(rows),length(cols))
#'       if(ivec[i]==jvec[i]){
#'         rows <- rows[row(tmat)[upper.tri(tmat)]]
#'         cols <- cols[col(tmat)[upper.tri(tmat)]]
#'       }else{
#'         rows <- rows[row(tmat)]
#'         cols <- cols[col(tmat)]
#'       }
#'       rows <- rows[nS!=0]
#'       cols <- cols[nS!=0]
#'       nS <- nS[nS!=0]
#'       cbind(rows,cols,nS)
#'     })
#'     gc()
#'   }
#'   return(tl)
#' }
#'
# fast_LD <- function(H,cummap,m,Ne,cutoff,nchunks){
 #  require(Matrix)
#'   require(matrixStats)
#'   nSNPs <- length(cummap)
#'   chunksize <-ceiling(nSNPs/nchunks)
#'   nSl <-fast_nS(H,cummap,m,Ne,cutoff,nchunks)
#'   rS <- Matrix(0,nrow = nSNPs,ncol = nSNPs,sparse = T)
#' for(i in 1:length(nSl)){
#'   cat(paste0(i,"\n"))
#'   tm <- values(nSl[[i]])
#'   rS[tm[tm[,"nS"]!=0,c("rows","cols")]] <-tm[tm[,"nS"]!=0,"nS"]
#'   rS[tm[tm[,"nS"]!=0,c("cols","rows")]] <-tm[tm[,"nS"]!=0,"nS"]
#' }
#'   diag(rS) <- colVars(H)
#'   nmsum <- sum(1/1:(2*m-1))
#'   theta <- (1/nmsum)/(2*m+1/nmsum)
#'   sighat <- (1-theta)^2*rS+0.5*theta*(1-0.5*theta)*Diagonal(nSNPs)
#'   cmat <- cov2cor(sighat)
#'   return(cmat)
#' }












#' Helper function for generating LD matrix
#'
#' @param Svec vector containing lower diagonal of reference panel covariance matrix
#' @param distvec vector containing lower diagonal of matrix of pairwise genetic distances
#' @param m integer sample size of genetic map (not of reference panel)
#' @param Ne numeric effective population size estimate for population
#' @param cutoff numeric cutoff for shrinkage
#' @param nSNPs integer number of SNPs
#' @param Sdiag vector containing diagonal of ref panel covariance matrix
#'
lddiag <- function(Svec,distvec,m,Ne,cutoff,nSNPs,Sdiag){
  nmsum <- sum((1/1:(2*m-1)))
  theta <- (1/nmsum)/(2*m+1/nmsum)

  rho <- 4*Ne*(distvec)/100
  shrinkage <- exp(-rho/(2*m))
  shrinkage[shrinkage<cutoff] <- 0
  nS <- shrinkage*Svec
  rS <- matrix(data = 0,nrow = nSNPs,ncol = nSNPs)
  rS[lower.tri(rS)] <- nS
  rS <- rS+t(rS)
  diag(rS)<- Sdiag
  sighat <- (1-theta)^2*rS+0.5*theta*(1-0.5*theta)*diag(nSNPs)
  return(sighat)
}


#'
#'
#' #' Helper function for generating LD matrix
#' #' Much slower than lddiag, but is as close as possible to the original MATLAB implementation
#' #' @param Svec vector containing lower diagonal of reference panel covariance matrix
#' slow.LD <- function(Hpanel,cummap,m,ne,cutoff){
#'   nmsum <- sum((1/1:(2*m-1)))
#'   theta <- (1/nmsum)/(2*m+1/nmsum)
#'   S <- covar(H)
#'   S[lower.tri(S)] <- 0
#'   nSNPs <- nrow(S)
#'   for(i in 1:nSNPs){
#'     j <- i+1
#'     while(j<=nSNPs){
#'       rho <- 4*Ne*(cummap[j]-cummap[i])/100
#'       shrinkage <- exp(-rho/(2*m))
#'       if(shrinkage<cutoff){
#'         shrinkage <- 0
#'       }
#'       S[i,j] <- shrinkage*S[i,j]
#'       j <- j+1
#'     }
#'   }
#'   td <- diag(S)
#'   S <- S+t(S)
#'   diag(S) <- td
#'   SigHat <-(1-theta)^2*S+0.5%*%theta%*%(1-0.5%*%theta)%*%diag(nSNPs)
#'   cmat <- cov2cor(SigHat)
#'   return(cmat)
#' }
#'
#'
#' Generates LD matrix from reference panel data
#'
#' @param ref_list a list with an element H, which is the matrix of genotypes, and a vector
#' cummap which is a vector with the cumulative map distance
#' @param m integer sample size of genetic map (not of reference panel)
#' @param Ne numeric effective population size estimate for population
#' @param cutoff numeric cutoff for shrinkage
gen_LD <- function(ref_list,m,Ne,cutoff){
  require(coop)
  H <- ref_list[["H"]]
  cummap <- ref_list[["cummap"]]
  distmat <- outer(cummap,cummap,'-')
  distvec <- distmat[lower.tri(distmat)]
  S <- covar(H)
  Svec <- S[lower.tri(S)]
  Sdiag <- diag(S)
  nSNPs <- nrow(S)
  ldmat <-lddiag(Svec,distvec,m,Ne,cutoff,nSNPs,Sdiag)
}
#' Generates LD matrix from reference panel data
#'
#' @param H, the matrix of genotypes,
#' @param cummap which is a vector with the cumulative map distance
#' @param m integer sample size of genetic map (not of reference panel)
#' @param Ne numeric effective population size estimate for population
#' @param cutoff numeric cutoff for shrinkage
#' @param chunksize size of chunks to use
arm_gen_LD <- function(H,cummap,m,Ne,cutoff,chunksize){
  require(future)
  plan(multiprocess)
  nSNPs <- ncol(H)
  nchunks <-ceiling(nSNPs/chunksize)
  resl <- list()
  for(i in 0:(nchunks-1)){
    for(j in i:(nchunks-1)){
      cat(paste0(i,"_",j,"\n"))
      resl[[paste0(i,"_",j)]] <-future({
        p_sparse_LD(cummap, H,Ne,m,cutoff,chunksize,i,j)
      })
    }
  }
  tm <- value(resl[[1]])
  fmat <-Reduce("+",lapply(resl,values))
  return(fmat)
}





