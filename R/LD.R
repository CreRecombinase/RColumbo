
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
  stopifnot(file.exists(haph5),file.exists(mapfile))
  mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
  legdf <- read_h5_df(haph5,"Legend")
  nSNPs <- nrow(legdf)
  if(length(rsids)==0){
  rslistvec <- legdf$rsid
  }else{
    rslistvec <-rsids
  }
  rslistvec <-intersect(rslistvec,legdf$rsid)
  rslistvec <- intersect(rslistvec,mapdf$rsid)
  stopifnot(length(rslistvec)>0)
  readind = which((legdf$rsid %in%rslistvec)&!duplicated(legdf$rsid))
  hapd <- read_haplotype_ind_h5(hap_h5file = haph5,indexes = readind)
  legdf <- slice(legdf,readind)
  mapdf <- filter(mapdf,rsid %in% rslistvec)
  colnames(hapd) <-legdf$rsid
  rhapd <- flipmat(hapd,as.integer(legdf$allele0<legdf$allele1))
  retlist <-list(H=rhapd,cummap=mapdf$cummap)
  return(retlist)
}



subset_rpanels <- function(h5file,mapfile){
  require(rhdf5)
  require(dplyr)
  mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
  legrsid <- h5vec(h5file,"Legend","rsid")
  return(intersect(mapdf$rsid,legrsid))
}

#'
#' #' Compute the intersection of genetic map,hdf5 file, and haplotype data for a particular chromosome
#' #' @param legendfile IMPUTE format legend file
#' #' @param hapfile IMPUTE format genotype file
#' #' @param mapfile file containing physical map info
#' #' @param h5file HDF5 formatted data with merged eQTL and GWAS data
subset_h5_ref <- function(haph5,mapfile,eqtlfile){
  require(readr)
  snpvec <- h5vec(eqtlfile,"SNP","rsid")
  snpdf <- data_frame(rsid=paste0("rs",snpvec))
  lsnpvec <- h5vec(haph5,"Legend","rsid")
  legenddf <- data_frame(rsid=lsnpvec)
  mapdat <- read.table(mapfile,sep=" ",header=F) %>% rename(rsid=V1,pos=V2,map=V3)
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


  # eqtlc <-group_by(eqtldf,Gene)%>% summarise(neqtl=n(),prode=sum(log(pvale),na.rm = T),prodg=sum(log(pvalg),na.rm=T))
  # eqtlc <- mutate(eqtlc,hash=as.integer(neqtl>=100))%>%arrange(desc(hash),prode)
  eqtldf <- filter(eqtldf,Gene==gene)

  legenddf <- read_delim(legendfile,delim = " ",col_names = T)
  #  hapmat <- data.matrix(read_delim(hapfile,delim = " ",col_names = F))
  mapdat <- read.table(mapfile,sep=" ",col.names=c("rsid","pos","dist"))
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
gen_LD <- function(ref_list,m=85,Ne=11490.672741,cutoff=1e-3){

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


intersect_df <- function(query_df,haph5,gwash5){
  gwas_rsid <- c(h5read(gwash5,"/GWAS/rsid"))
  query_df <- filter(query_df,rsid %in% gwas_rsid)
  leg_df <- read_h5_df(haph5,"Legend") %>% mutate(haplotype_ind=1:length(rsid))
  map_df <- read_h5_df(haph5,"Map")
  nsnpdf <- inner_join(map_df,leg_df) %>% inner_join(select(query_df,-rsid,-doFlip)) %>% filter(!duplicated(rsid))
  return(nsnpdf)
}


# fast_dense_LD <- function(h5file,offset,chunksize,m=85,Ne=11490.672741,cutoff=1e-3){
#
# }

 #Return Dense LD matrix from character vector of rsid's
comp_dense_LD <- function(indvec,cummap,rsid,haph5,m=85,Ne=11490.672741,cutoff=1e-3){
  require(Matrix)
  require(dplyr)
  stopifnot(file.exists(haph5),length(cummap)==length(indvec),length(indvec)==length(rsid))
  fmat <- gen_dense_LD(haph5,indvec,cummap,m,Ne,cutoff,0,0,length(rsid))
  fmat <- fmat+t(fmat)
  diag(fmat)<- 1
  rownames(fmat) <-rsid
  colnames(fmat) <- rsid
  return(fmat)
}



gen_LD_rsid <- function(rsidvec,haph5f,m=85,Ne=11490.672714,cutoff=1e-3){
  require(rhdf5)
  rsidfilt <- function(x){x %in% rsidvec}
  full_rsid <- h5read(haph5f,"/Legend/rsid")
  which_rsid <- which(full_rsid %in% rsidvec)
  pos_df <- read_h5_df_filter(haph5f,"Legend","rsid",filter_fun = rsidfilt)
  mapdf <- read_h5_df_filter(haph5f,"Map","rsid",filter_fun=rsidfilt)
  leg_map <- inner_join(mapdf,pos_df)
  cummap <- leg_map$cummap
  stopifnot(nrow(leg_map)==length(which_rsid),length(which_rsid)==length(rsidvec))
  fmat <- comp_dense_LD(which_rsid,cummap,rsidvec,haph5f,m,Ne,cutoff)
  return(fmat)
}







#' Generates LD matrix from reference panel data
#'
#' @param rsids vec
#' @param cummap which is a vector with the cumulative map distance
#' @param m integer sample size of genetic map (not of reference panel)
#' @param Ne numeric effective population size estimate for population
#' @param cutoff numeric cutoff for shrinkage
#' @param chunksize size of chunks to use
torque_arm_gen_LD <- function(rslistvec=character(0),eqtlfile=tempfile(),haph5,mapfile,m=85,Ne=11490.672741,cutoff=1e-3,chunksize,i,j,result_dir,chrom){
  require(Matrix)
  require(dplyr)
  stopifnot(file.exists(mapfile),
            file.exists(haph5),
            (file.exists(eqtlfile)||length(rslistvec)>0))
  doSkip<- FALSE
  if(j>0){
    tj <- j-1
    prev_outfile <-file.path(result_dir,paste0("chr",chrom),paste0("chr",chrom,"_",i,"_",tj,"_",chunksize,".RDS"))
    if(file.exists(prev_outfile)){
      nmat <- readRDS(prev_outfile)
      if(sum(nmat)==0){
        doSkip <- TRUE
      }
    }
  }

  mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
  legdf <- read_h5_df(haph5,"Legend")
  if(length(rslistvec)==0){
    rslistvec <- subset_h5_ref(haph5 = haph5,mapfile = mapfile,eqtlfile = eqtlfile)
    stopifnot(length(rslistvec)>0)
    rslistvec <-intersect(rslistvec,legdf$rsid)
    rslistvec <- intersect(rslistvec,mapdf$rsid)
    stopifnot(length(rslistvec)>0)
  }
  readind = which((legdf$rsid %in%rslistvec)&!duplicated(legdf$rsid))
  legdf <- slice(legdf,readind)
  mapdf <- filter(mapdf,rsid %in% rslistvec)
  cummap <- mapdf$cummap
  nSNPs <- length(cummap)
  nchunks <-ceiling(nSNPs/chunksize)
  rm(mapdf,legdf,rslistvec)
  if(!doSkip){
    fmat <- flip_hap_LD(haph5,readind,cummap,m,Ne,cutoff,i,j,chunksize)
  }else{
    fmat <- sparseMatrix(i = 1,j=1,x = 0,dims=c(nSNPs,nSNPs))
  }
  outfile <-file.path(result_dir,paste0("chr",chrom),paste0("chr",chrom,"_",i,"_",j,"_",chunksize,".RDS"))
  saveRDS(fmat,outfile)
  return(outfile)
}


fhLD <- function(haph5,readind,doFlip,cummap,m,Ne,cutoff,i,j,chunksize){

  nmsum <- sum((1/1:(2*m-1)))
  theta <- (1/nmsum)/(2*m+1/nmsum)
  readind = which((totlegdf$rsid %in%rslistvec)&!duplicated(totlegdf$rsid))
  nSNPs <- length(cummap)
  nchunks <-ceiling(nSNPs/chunksize)
  mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
  legdf <- read_h5_df(haph5,"Legend")
  readind = which((legdf$rsid %in%rslistvec)&!duplicated(legdf$rsid))
  legdf <- slice(legdf,readind)

  doFlip <- as.integer(legdf$allele0<legdf$allele1)
  istart=i*chunksize+1
  jstart=j*chunksize+1
  istop=min((i+1)*chunksize,nSNPs)
  jstop=min((j+1)*chunksize,nSNPs)

  mapa= cummap[istart:istop]
  mapb= cummap[jstart:jstop]

  hmata <- flip_hap(haph5,readind,doFlip,i,chunksize,nSNPs)
  hmatb <- flip_hap(haph5,readind,doFlip,j,chunksize,nSNPs)

  distmat <-outer(mapb,mapa,"-")
  if(i==j){
    distmat[upper.tri(distmat)]<- 0
    distmat <- distmat+t(distmat)
    distmat[lower.tri(distmat)] <- 0
  }

  S <- cov(hmata,hmatb)
  if(i==j){
    S[lower.tri(S)]<-0
  }
  distmat=4*Ne*distmat/100
  distmat=exp(-distmat/(2*m))
  distmat[distmat<cutoff] <- 0
  distmat=distmat*S
  if(i==j){
    diag(distmat) <- apply(hmata,2,var)
    distmat=(1-theta)*(1-theta)*distmat+0.5*theta*(1-0.5*theta)*diag(nrow(distmat))
    vars=diag(distmat)
    vars =1/sqrt(vars);
    distmat <- sweep(distmat,1,vars,"*")
    distmat <- sweep(distmat,2,vars,"*")
    diag(distmat) <- 1
  }
  return(distmat)
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

subset_all_hap <- function(query_df,haph5files,gwash5){
  require(dplyr)
  gwas_rsid <- c(h5read(gwash5,"/GWAS/rsid"))
  query_df <- filter(query_df,rsid %in% gwas_rsid)
  snpl <- list()
  for(i in 1:22){
    cat(paste0(i,"\n"))
    thf <- read_h5_df(haph5files[i],"Legend") %>% mutate(chrom=i,hsnpid=1:length(rsid))
    tmf <- read_h5_df(haph5files[i],"Map") %>% mutate(chrom=i)
    snpl[[i]] <-inner_join(thf,tmf)
  }
  asnpdf <- bind_rows(snpl)
  nsnpdf <- inner_join(asnpdf,select(query_df,-rsid,-doFlip)) %>% filter(!duplicated(rsid))
  return(nsnpdf)
}


sslab.em <- function(p,bh,si){
  pi <- p[1]
  tau <- p[2]
#  uzin <- pi*dnorm(bh,mean=0,sd=sqrt(tau^2+si^2))
  #uizd <- uzin+(1-pi)*dnorm(bh,mean=0,sd=si,log=F)
  #uiz <- uzin/uizd
  uiz <- 1/(1+((1-pi)*dnorm(bh,mean=0,sd=si,log=F))/(pi*dnorm(bh,mean=0,sd=sqrt(tau^2+si^2))))

  # cat(paste0("pi:",pi,"\n"))
  # cat(paste0("tau:",tau,"\n"))
  nt <- sqrt(sum(uiz*bh^2)/sum(uiz)-sum(si^2*uiz)/sum(uiz))
  np <- mean(uiz)
  return(c(np,nt))
}

sslab.lik <- function(p,bh,si){
  pi <- p[1]
  tau <- p[2]
  uzin <- pi*dnorm(bh,mean=0,sd=sqrt(tau^2+si^2))
  uizd <- uzin+(1-pi)*dnorm(bh,mean=0,sd=si,log=F)
  uiz <- uzin/uizd
  return(-sum(uiz*dnorm(bh,mean=0,sd=sqrt(tau^2+si^2))+(1-uiz)*dnorm(bh,mean=0,sd=si)))
}


independent_betahat <- function(h5file,LDcutoff,chunksize=50000){
  require(rhdf5)
  require(dplyr)
  require(tidyr)
  require(BBmisc)
  require(SQUAREM)
  mfmat <- chunk_LD(h5file = h5file,offset = 0,chunksize = chunksize,m=85,Ne=11490.672741,cutoff=1e-3)


}




estimate_pi_tau <- function(rawh5,sub_snpleg,subleg,outh5,chunksize=100000){
  require(rhdf5)
  require(dplyr)
  require(tidyr)
  require(BBmisc)
  require(SQUAREM)
  stopifnot(nrow(subleg)==1)
  chunkleg <- chunk_df(sub_snpleg,n.chunks = 2)
  expdat <- read_dmat_ind_h5(rawh5,"EXPdata","orthoexpression",c(subleg$exp_id))
  expsd <- apply(expdat,2,sd)
  n <- nrow(expdat)
  gp <- c(runif(1),exp(runif(1)))
  betavl <- list()
  serrvl <- list()
  for(i in 1:length(chunkleg)){
    cat(paste0("chunk: ",i," of ",length(chunkleg),"\n"))
    tchunk <- chunkleg[[i]]
    rawindvec <- tchunk$snp_ind
    osnpdat <- read_dmat_chunk_ind(rawh5,"SNPdata","orthogenotype",rawindvec)
    snpsd <- c(colssd(osnpdat))
    rmat <- rMatrix(osnpdat,expdat)
    betav <-c((rmat*expsd)/snpsd)
    betavl[[i]] <- betav
    serrv <- c(betav/(sqrt(n-2)*rmat/sqrt(1-rmat^2)))
    serrvl[[i]] <- serrv
  }
  bh <- unlist(betavl)
  si <- unlist(serrvl)
  emr <- try(squarem(par=gp,bh=bh,si=si,fixptfn=sslab_em,control=list(trace=T)))
  if(inherits(emr,"try-error")){
    j <- 0
    while(inherits(emr,"try-error")){
      cat(paste0(j,"\n"))
      tgp <- c(runif(1),exp(runif(1)))
      emr <- try(squarem(par = tgp,y=y,fixptfn = sslab.em,objfn = sslab.lik,control = list(trace=T)))
      j <- j+1
    }
  }
  fpi <- emr$par[1]
  ftau <- emr$par[2]
  data <- data_frame(bh=bh,si=si)
  data <- mutate(data,mu=c(pMu(bh,si,pi=fpi,tau=ftau)))
  data <- mutate(data,rsidi=sub_snpleg$rsidi)
  sub_snpleg <- inner_join(sub_snpleg,data)
  sub_snpleg <- mutate(sub_snpleg,isCis=chrom==subleg$chrom&((abs(pos-subleg$TSStart)<1e6|abs(pos-subleg$TSStop)<1e6)|(pos<subleg$TSStop&pos>subleg$TSStart)))



  ggplot(sub_snpleg)+geom_histogram(aes(x=mu,..density..),binwidth = 0.01)+facet_wrap(~isCis)



  return(emr)
}






stat_extract <- function(rawh5,haph5,gwash5,outh5,chromosome,chunksize,cis_pcutoff=0.01,trans_pcutoff,cis_LDcutoff,trans_LDcutoff,cisdist_cutoff=1e6,append=F){
  require(rhdf5)
  require(dplyr)
  require(tidyr)
  require(BBmisc)
  # eqtlh5 <- "~/Desktop/eQTL/Snake/IBD_WholeBlood_eQTL.h5"
  # rawh5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_raw_data.h5"
  # haph5 <- "/home/nwknoblauch/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno_hap.h5"
  # gwash5 <- "~/Desktop/eQTL/Snake/IBD.h5"
  cat(paste0("Preparing data for chromosome:",chromosome,"\n"))
  snpleg <- read_h5_df(rawh5,"SNPinfo") %>% mutate(snp_ind=1:n())
  expleg <-read_h5_df(rawh5,"EXPinfo") %>% mutate(exp_ind=1:n())
  bexpdat <- read_dmat_ind_h5(rawh5,"EXPdata","orthoexpression",c(1:23973))
  cis_tcutoff <- abs(qt(cis_pcutoff,df = nrow(bexpdat)-2,lower.tail = F))
  trans_tcutoff <- abs(qt(trans_pcutoff,df=nrow(bexpdat)-2,lower.tail=F))
  ch19leg <- filter(snpleg,chrom==chromosome)
  inter_df <- intersect_df(ch19leg,haph5,gwash5 = gwash5)
  dfl <- chunk_df(inter_df,chunk.size = chunksize)
  chunkind <- chunk(1:(nrow(inter_df)),chunk.size = chunksize)
  if(file.exists(outh5)&(!append)){
    file.remove(outh5)
  }
  for(i in 1:length(dfl)){
    cat(paste0(i," of ",length(dfl),"\n"))
    query_df <- dfl[[i]]
    snpi <- query_df$snp_ind
    query_df <- mutate(query_df,rsid=as.integer(gsub("rs","",rsid)))
    rsids <- query_df$rsid
    expids <- expleg$fgeneid

    fmat <- comp_dense_LD(query_df$haplotype_ind,query_df$cummap,query_df$rsid,haph5,m=85,Ne=11490.672741,cutoff=1e-3)
    osnpdat <- read_dmat_ind_h5(rawh5,"SNPdata","orthogenotype",snpi)

    # betas <- betaMatrix(osnpdat,bexpdat)
    rmat <- rMatrix(osnpdat,bexpdat)

    cis_eqtl <- extract_stats(Genotype = osnpdat,snpanno = query_df,Expression =  bexpdat,
                              expanno = expleg,LDmat = fmat,
                              rmat=rmat,tcutoff = cis_tcutoff,
                              LDcutoff = cis_LDcutoff,
                              cisdist = 1e6,
                              display_progress = T,doCis = T)
    trans_eqtl<- extract_stats(Genotype = osnpdat,snpanno = query_df,Expression =  bexpdat,
                               expanno = expleg,LDmat = fmat,
                               rmat=rmat,tcutoff = trans_tcutoff,
                               LDcutoff = trans_LDcutoff,
                               cisdist = 1e6,
                               display_progress = T,doCis = F)
    write_Rnumeric_h5(outh5,"cis_eQTL","thetahat",data = cis_eqtl$theta,deflate_level = 4)
    write_Rnumeric_h5(outh5,"cis_eQTL","serr",data = cis_eqtl$serr,deflate_level = 4)
    write_Rint_h5(outh5,"cis_eQTL","rsid",data = cis_eqtl$rsid,deflate_level = 4)
    write_Rint_h5(outh5,"cis_eQTL","fgeneid",data = cis_eqtl$fgeneid,deflate_level = 4)

    write_Rnumeric_h5(outh5,"trans_eQTL","thetahat",data = trans_eqtl$theta,deflate_level = 4)
    write_Rnumeric_h5(outh5,"trans_eQTL","serr",data = trans_eqtl$serr,deflate_level = 4)
    write_Rint_h5(outh5,"trans_eQTL","rsid",data = trans_eqtl$rsid,deflate_level = 4)
    write_Rint_h5(outh5,"trans_eQTL","fgeneid",data = trans_eqtl$fgeneid,deflate_level = 4)
  }
  cat("Done!\n")
}

chunk_eQTL <- function(exph5,snph5,outh5,snpinter=NULL,expinter=NULL,cisdist_cutoff=1e6){
  require(dplyr)
  require(tidyr)
  require(BBmisc)
  require(h5)
  cat("Starting to read snph5:",snph5)
  snpleg <-read_df_h5(snph5,"SNPinfo",filtervec=snpinter)
  cat("Starting to read exph5:",exph5)
  expleg <- read_df_h5(exph5,"EXPinfo",filtervec=expinter)
  expdat <- read_dmat_chunk_ind(exph5,"EXPdata","expression",expinter)
  snpdat <- read_dmat_chunk_ind(snph5,"SNPdata","genotype",snpinter)
#  feqtl <- really_fast_eQTL(Genotype = snpdat,snpanno = snpleg,Expression = expdat,expanno = expleg)
  eqtl <- fast_eQTL(Genotype=snpdat,snpanno=snpleg,Expression=expdat,
                    expanno=expleg,cis_tcutoff = 0,trans_tcutoff = 0,cisdist = 1e6,doTrans = T,doCis = T)
  cis_eqtl <-filter(eqtl,cistrans==1) %>% select(-cistrans) %>% mutate(tstat=theta/serr)
  trans_eqtl <- filter(eqtl,cistrans==0) %>% select(-cistrans) %>% mutate(tstat=theta/serr)
  if(nrow(cis_eqtl)>0){
    write_h5_df(cis_eqtl,"cis_eQTL",outfile = outh5)
  }
  if(nrow(trans_eqtl)>0){
    write_h5_df(trans_eqtl,"trans_eQTL",outfile = outh5)
  }
  gc()
  cat("Done!\n")
  return(nrow(eqtl))
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
  eqtl <- fastest_eQTL(genotypef=snph5,snpinter=snpinter,
                       expressionf=exph5,expinter=expinter)
  write_2dmat_h5(h5f = outh5,groupn = "eQTL",datan = "beta_mat",chunksize = c(length(snpinter)/2,length(expinter)/2),deflate_level = 4,data = eqtl[1,,])
  write_2dmat_h5(h5f = outh5,groupn = "eQTL",datan = "t_mat",chunksize = c(length(snpinter)/2,length(expinter)/2),deflate_level = 4,data = eqtl[2,,])

  cat("Done!\n")
  return(dim(eqtl))
}


run_eqtl<- function(rawh5,outh5,chromosome,chunksize,cis_pcutoff=0.01,trans_pcutoff=1e-3,cisdist_cutoff=1e6,append=F,useortho=F){
  require(dplyr)
  require(tidyr)
  require(BBmisc)
  # eqtlh5 <- "~/Desktop/eQTL/Snake/IBD_WholeBlood_eQTL.h5"
  # rawh5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_raw_data.h5"
  # haph5 <- "/home/nwknoblauch/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno_hap.h5"
  # gwash5 <- "~/Desktop/eQTL/Snake/IBD.h5"
  cat(paste0("Preparing data for chromosome:",chromosome,"\n"))
  snpleg <- read_h5_df(rawh5,"SNPinfo") %>% mutate(snp_ind=1:n())
  expleg <-read_h5_df(rawh5,"EXPinfo") %>% mutate(exp_ind=1:n())
  # expleg <- as_data_frame(t(h5read(rawh5,"EXPinfo/annomat"))) %>% rename(fgeneid=V1,sgeneid=V2,chrom=V3,TSStart=V4,TSStop=V5) %>% mutate(exp_id=0:(length(fgeneid)-1))
  # expleg <- rename(expleg,start=TSStart,end=TSStop)
  if(useortho){
    bexpdat <- read_dmat_ind_h5(rawh5,"EXPdata","orthoexpression",c(1:get_rownum_h5(rawh5,"EXPdata","orthoexpression")))
  }else{
    bexpdat <- read_dmat_ind_h5(rawh5,"EXPdata","expression",c(1:get_rownum_h5(rawh5,"EXPdata","expression")))
    bexpdat <- scale(bexpdat,center = T,scale=T)
  }




  cis_tcutoff <- abs(qt(cis_pcutoff,df = nrow(bexpdat)-2,lower.tail = F))

  trans_tcutoff <- abs(qt(trans_pcutoff,df=nrow(bexpdat)-2,lower.tail=F))
  chleg <- filter(snpleg,chrom==chromosome)
  dfl <- chunk_df(chleg,chunk.size = chunksize)
  chunkind <- chunk(1:(nrow(chleg)),chunk.size = chunksize)
  if(file.exists(outh5)&(!append)){
    file.remove(outh5)
  }
  for(i in 1:length(dfl)){
    cat(paste0(i," of ",length(dfl),"\n"))
    query_df <- dfl[[i]]
    snpi <- query_df$snp_ind
    # fmat <- comp_dense_LD(query_df$haplotype_ind,query_df$cummap,query_df$rsid,haph5,m=85,Ne=11490.672741,cutoff=1e-3)
    # tsnp <- osnpdat[,13]
    # texp <- bexpdat[,1]
    # summary(lm(texp~tsnp+0))
    if(useortho){
      osnpdat <- read_dmat_ind_h5(rawh5,"SNPdata","orthogenotype",snpi)
    }
    else{
      osnpdat <- read_dmat_ind_h5(rawh5,"SNPdata","genotype",snpi)
    }
    # betas <- betaMatrix(osnpdat,bexpdat)
    eqtl <- fast_eQTL(Genotype=osnpdat,snpanno=query_df,Expression=bexpdat,
                      expanno=expleg,cis_tcutoff = cis_tcutoff,trans_tcutoff = trans_tcutoff,cisdist = 1e6,doTrans = T,doCis = T)
    cis_eqtl <-filter(eqtl,cistrans==1) %>% select(-cistrans) %>% mutate(tstat=theta/serr)
    trans_eqtl <- filter(eqtl,cistrans==0) %>% select(-cistrans)
    if(nrow(cis_eqtl)>0){
      write_h5_df(cis_eqtl,"cis_eQTL",outfile = outh5)
    }
    if(nrow(trans_eqtl)>0){
      write_h5_df(trans_eqtl,"trans_eQTL",outfile = outh5)
    }
    gc()
  }
  cat("Done!\n")
}






