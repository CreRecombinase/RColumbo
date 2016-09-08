#' Reads an HDF5 file into dataframes according to the size of each vector, and the output of
#' @param lsdf dataframe (as output from h5ls)
#' @param h5file location of h5file
h5df <- function(gdf,h5file){
  stopifnot(all(gdf$dim==gdf$dim[1]))
  lsn <- paste0(gdf$group,"/",gdf$name)
  datl <- lapply(lsn,function(x,filename){
    tn <- strsplit(x,"/",fixed=T)[[1]][3]
    tdat <- c(h5read(file = filename,name = x))
    tdf = data_frame(x=tdat)
    colnames(tdf) <- tn
    return(tdf)
  },filename=h5file)
  hdf <- bind_cols(datl)
  return(hdf)
}
h5vec <-function(h5file,groupname,dataname){
  tn <-paste0("/",groupname,"/",dataname)
  retvec <- c(h5read(h5file,tn))
  return(retvec)
}

#' Turns matrix of diagonals in to sparse banded matrix
#'
#' @param bmat matrix of diagonals
#' @param bwd band width
from_band <- function(bmat){
  tmat <- bandSparse(n = ncol(bmat),m = ncol(bmat),k=-c(0:(nrow(bmat)-1)),diagonals = t(bmat),symmetric = T)
  return(tmat)
}

readMATLAB <- function(matfile,group="BR"){
  require(rhdf5)
  matres <- from_band(h5read(matfile,name=group))
}



matLines <- function(gzc,nr,nc){
  tL <-readLines(gzc,n = nr)
  tLL <- matrix(as.integer(unlist(strsplit(tL,split = " ",fixed = T,useBytes = T))),
                nrow = nr,ncol = nc,byrow = T)
  return(tLL)
}





#' Write haplotype data to HDF5 data
#' @param haph5 h5file (with legend data already written)
#' @param chunksize chunk size for compression
write_doFlip <- function(haph5,chunksize=100000){
  legdat <- read_h5_df(haph5,"Legend")
  legdat$doFlip <- as.integer(legdat$allele0<legdat$allele1)
  if(sum(colnames(legdf)=="doFlip")==0){
    h5createDataset(haph5,"/Legend/doFlip",
                    dims=c(nrow(legdat)),
                    storage.mode="integer",
                    chunk=chunksize,level=2)
    chunkind <- chunk(1:nSNP,chunk.size = chunksize)
    for(i in 1:length(chunkind)){
      cat(paste0(i,"\n"))
      tldat <- slice(legdat,chunkind[[i]])
      h5write(tldat$doFlip,file=haph5,name="/Legend/doFlip",index=list(chunkind[[i]]))
    }
  }
}


#'Write covariate data to HDF5 data
#'@param covarf covariate file
#'@param h5file output H5file
#'@param chunksize (this is not really necessary)
write_covar_h5 <- function(covarf,h5file,chunksize=1,deflate_level=9){
  require(rhdf5)
  library(dplyr)
  covardat <- read.table(covarf,header=T,stringsAsFactors = F)
  matdat <-t(data.matrix(select(covardat,-ID)))
  ncovar <- ncol(matdat)
  nid <- nrow(matdat)
  tmatdat <-matdat[,1:ceiling(ncovar/2)]
  write_dmatrix_h5(h5file,"Covardat", "covariates",ncovar, nid, tmatdat,deflate_level = deflate_level)
  tmatdat <-matdat[,(ceiling(ncovar/2)+1):ncol(matdat)]
  write_dmatrix_h5(h5file,"Covardat", "covariates",ncovar, nid, tmatdat,deflate_level = deflate_level)
  h5createGroup(file = h5file,group = "Covarinfo")
  h5createDataset(h5file,"/Covarinfo/id",
                  dims=c(nrow(covardat)),
                  storage.mode="character",size=max(nchar(covardat$ID))+1,
                  chunk=chunksize,level=deflate_level)
  h5write(covardat$ID,file=h5file,name="/Covarinfo/id")

}


#' Write haplotype data to HDF5 data
#' @param haplotype matrix (one SNP per column, one individual per row)
#' @param legdat legend dataframe
#' @param haph5 output HDF5 filename
#' @param chunksize chunk size for compression
write_leg_h5 <- function(legendfile,haph5,chunksize=100000){
  require(rhdf5)
  library(BBmisc)
  library(dplyr)
  legdat <- read.table(legendfile,header=T,sep=" ",stringsAsFactors=F)
  legdat$doFlip <- as.integer(legdat$a0<legdat$a1)
  nind <- length(scan(hapfile,what=integer(),sep=" ",nlines=1))
  nSNP <- nrow(legdat)
  chunkind <- chunk(1:nSNP,chunk.size = chunksize)
  maxchar <- max(nchar(legdat$id))+1
  if(!file.exists(haph5)){
    h5createFile(haph5)
  }

  h5createGroup(file = haph5,group = "Legend")
  h5createDataset(haph5,"/Legend/doFlip",
                  dims=c(nrow(legdat)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Legend/rsidi",
                  dims=c(nrow(legdat)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Legend/rsid",
                  dims=c(nrow(legdat)),
                  storage.mode="character",
                  size=maxchar+1,
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Legend/pos",
                  dims=c(nrow(legdat)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Legend/allele0",
                  dims=c(nrow(legdat)),
                  storage.mode="character",
                  size=max(nchar(legdat$a0))+1,
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Legend/allele1",
                  dims=c(nrow(legdat)),
                  storage.mode="character",
                  size=max(nchar(legdat$a1))+1,
                  chunk=chunksize,level=2)
  for(i in 1:length(chunkind)){
    cat(paste0(i,"\n"))
    tldat <- slice(legdat,chunkind[[i]])
    tldat$rsidi <- as.integer(gsub("rs","",tldat$id))
    tldat$rsidi <- ifelse(is.na(mapdf$rsidi),0,tldat$rsidi)
    h5write(tldat$rsidi,file=haph5,name="/Legend/rsidi",index=list(chunkind[[i]]))
    h5write(tldat$doFlip,file=haph5,name="/Legend/doFlip",index=list(chunkind[[i]]))
    h5write(tldat$id,file=haph5,name="/Legend/rsid",index=list(chunkind[[i]]))
    h5write(tldat$position,file=haph5,name="/Legend/pos",index=list(chunkind[[i]]))
    h5write(tldat$a0,file=haph5,name="/Legend/allele0",index=list(chunkind[[i]]))
    h5write(tldat$a1,file=haph5,name="/Legend/allele1",index=list(chunkind[[i]]))
  }
    H5close()
}

rewrite_h5 <- function(haph5,chunksize=100000){
  hrsidi <- as.integer(gsub("rs","",c(h5read(haph5,"/Legend/rsid"))))
  hrsidi <- ifelse(is.na(hrsidi),0,hrsidi)
  mrsidi <- as.integer(gsub("rs","",c(h5read(haph5,"/Map/rsid"))))
  mrsidi <- ifelse(is.na(mrsidi),0,mrsidi)
  H5close()
  h5createDataset(haph5,"/Legend/rsidi",
                  dims=c(length(hrsidi)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Map/rsidi",
                  dims=c(length(mrsidi)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5write(hrsidi,file=haph5,name="/Legend/rsidi")
  h5write(mrsidi,file=haph5,name="/Map/rsidi")
  H5close()
}


write_map_h5 <- function(mapfile,haph5,chunksize=125000){
  require(rhdf5)
  library(BBmisc)
  require(dplyr)
  mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
  nSNP <- nrow(mapdf)
  chunkind <- chunk(1:nSNP,chunk.size = chunksize)
  maxchar <- max(nchar(mapdf$rsid))
  if(!file.exists(haph5)){
    h5createFile(haph5)
  }
  h5createGroup(file = haph5,group = "Map")
  h5createDataset(haph5,"/Map/rsid",
                  dims=c(nrow(mapdf)),
                  storage.mode="character",
                  size=maxchar,
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Map/pos",
                  dims=c(nrow(mapdf)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Map/rsidi",
                  dims=c(nrow(mapdf)),
                  storage.mode="integer",
                  chunk=chunksize,level=2)
  h5createDataset(haph5,"/Map/cummap",
                  dims=c(nrow(mapdf)),
                  storage.mode="double",
                  chunk=chunksize,level=2)
  mapdf$rsidi <- as.integer(gsub("rs","",mapdf$rsid))
  mapdf$rsidi <- ifelse(is.na(mapdf$rsidi),0,mapdf$rsidi)
h5write(mapdf$rsidi,file=haph5,name="/Map/rsidi")
    h5write(mapdf$rsid,file=haph5,name="/Map/rsid")
    h5write(mapdf$pos,file=haph5,name="/Map/pos")
    h5write(mapdf$cummap,file=haph5,name="/Map/cummap")
}


read_h5_df_filter <- function(h5file,groupname,colname,filter_fun){
  require(rhdf5)
  idcol <- c(h5read(h5file,paste0("/",groupname,"/",colname)))
  indexes <- which(filter_fun(idcol))
  cols <- h5ls(h5file) %>% filter(group==paste0("/",groupname)) %>% select(name)
  cols <- cols$name
  retdf <- bind_cols(lapply(cols,function(x){
    td <- data_frame(c(h5read(h5file,paste0("/",groupname,"/",x),index=list(indexes))))
    colnames(td) <- x
    return(td)
  }))
  return(retdf)
}


write_h5_df <- function(df,group,outfile,deflate_level=4){
  dataname <- colnames(df)
  for(i in 1:length(dataname)){
    td <- df[[dataname[i]]]
    if(typeof(td)=="integer"){
      write_Rint_h5(h5file = outfile,groupname = group,dataname = dataname[i],data = td,deflate_level = deflate_level)
    }
    else{
      if(typeof(td)=="double"){
        write_Rnumeric_h5(h5file = outfile,groupname = group,dataname = dataname[i],data = td,deflate_level = deflate_level)
      }
      else{
        stop(paste0("type for ",dataname[i]," is :",typeof(td)," no matching method to write to HDF5"))
      }
    }
  }
}

read_h5_df <- function(h5file,groupname){
  require(rhdf5)
  require(dplyr)
  cols <- h5ls(h5file) %>% filter(group==paste0("/",groupname)) %>% select(name)
  cols <- cols$name
  retdf <- bind_cols(lapply(cols,function(x){
    td <- data_frame(c(h5read(h5file,paste0("/",groupname,"/",x))))
    colnames(td) <- x
    return(td)
  }))
  return(retdf)
}


#' Write haplotype data to HDF5 data
#' @param haplotype matrix (one SNP per column, one individual per row)
#' @param h5file file output file
#' @param chunksize
read_hap_h5 <- function(haph5,rslistvec=character(0),poslist=integer(0)){
  stopifnot(length(rslistvec)>0||length(poslist)>0)
  legdf <- read_h5_df(haph5,"Legend")
  nSNPs <- nrow(legdf)
  readind <- which((legdf$rsid %in% rslistvec)&!duplicated(legdf$rsid))
  legdf <- slice(legdf,readind)
  doFlip <- as.integer(legdf$allele0<legdf$allele1)
  return(flip_hap(haph5,readind,doFlip,0,length(readind),nSNPs))
}





#'Returns list of dataframes,one for each dimension of the data
#' @param h5file location of h5file
h5dfl <- function(h5file){
  require(rhdf5)
  require(dplyr)
  lsdf <- h5ls(h5file) %>% filter(otype!="H5I_GROUP")
  lsdl <- split(lsdf,lsdf$group)
  retl <- lapply(lsdl,h5df,h5file=h5file)
  names(retl) <- substring(names(retl),2)
  return(retl)
}

#' Merge VEGAS gene-based p-values with gencode file (to get gene symbol and ENSGID)
#' @param vegasf path to vegas results
#' @param gmf path to gencode file
vegas.merge <- function(vegasf,gmf){
 # eg <- rename(eg,thetahat=beta,eqtltstat=tstat,eqtlp=pval)
  gmd <- read.table(gmf,header=F,sep="\t",skip=5,quote='',stringsAsFactors = F)
  gmd <- mutate(gmd,gname=gsub(".+gene_name \"(.+)\"; transcript_type.+","\\1",V9))
  gmd <- mutate(gmd,ensid=gsub("gene_id \"(ENSG[0-9]+).[0-9]+\"; transcript_id.+","\\1",V9))
  gmd <- select(gmd,gname,ensid)
  gmd <- rename(gmd,Gene=gname)
  gmd <- distinct(gmd,Gene,ensid)
  vegd <- read.table(vegasf,header=T,stringsAsFactors = F,sep=" ")
  vegd <- arrange(vegd,Pvalue) %>% select(Gene,vegp=Pvalue,veg_snp_p=SNP.pvalue,veg_snps=nSNPs)

  gmd <- inner_join(gmd,vegd)
  return(gmd)
}

#' Merge eQTL results with gencode file to get gene symbols for every eQTL
#' @param eqtldf eQTL dataframe
#' @param gmf path to gencode file
gencode_eqtl <- function(eqtldf,gmf){
  # gmf <- "/media/nwknoblauch/Data/GTEx/gencode.v19.genes.patched_contigs.gtf.gz"
  gencode_df <- read.table(gmf,header=F,sep="\t",skip=5,quote='',stringsAsFactors = F)
  gencode_df <- mutate(gencode_df,gname=gsub(".+gene_name \"(.+)\"; transcript_type.+","\\1",V9))
  gencode_df <- mutate(gencode_df,ensid=gsub("gene_id \"(ENSG[0-9]+).[0-9]+\"; transcript_id.+","\\1",V9))
  gencode_df <- select(gencode_df,Gene=gname,ensid) %>% distinct(Gene,ensid)
  eqtldf <- mutate(eqtldf,ensid=paste0("ENSG",sprintf("%011d",fgeneid)))
  eqtldf <- inner_join(eqtldf,gencode_df,by="ensid")
  return(eqtldf)
}

sparse_cov2cor <-function(rdsfile,chunksize,variances,nSNPs){
  cat(rdsfile)
  i <- as.integer(gsub(paste0(".+chr.+_([0-9]+)_[0-9]+_",chunksize,".RDS"),"\\1",rdsfile))
  j <- as.integer(gsub(paste0(".+chr.+_([0-9]+)_([0-9]+)_",chunksize,".RDS"),"\\2",rdsfile))
  istart=(i*chunksize)+1;
  istop=(min((i+1)*chunksize-1,nSNPs-1))+1
  jstart=(i*chunksize)+1;
  jstop=(min((i+1)*chunksize-1,nSNPs-1))+1

  ivar <- 1/sqrt(variances[istart:istop])
  jvar <- 1/sqrt(variances[jstart:jstop])
  tdat <- readRDS(rdsfile)
  stdat <-as.matrix(tdat[istart:istop,jstart:jstop])
  if(i==j){
    tv <- diag(stdat)
    stdat <- stdat+t(stdat)
    diag(stdat) <-tv
    Is <-sqrt(1/diag(stdat))
    r <- stdat
    r[]
    r <- sweep(sweep(stdat,1,ivar,"*"),2,jvar,"*")
  }
  r <- sweep(sweep(stdat,1,ivar,"*"),2,jvar,"*")
  tdat[istart:istop,jstart:jstop] <-r
  return(tdat)
}




# mydat <- matrix(runif(30),6,5)
# V <- cov(mydat)
# Is <- sqrt(1/diag(V)) # diag( 1/sigma_i )
# r <- V # keep dimnames
# r[] <- Is * V * rep(Is, each = p)
# ##	== D %*% V %*% D  where D = diag(Is)
# r[cbind(1:p,1:p)] <- 1 # exact in diagonal
# r
# tcov2cor <- function(V)
# {
#   ## Purpose: Covariance matrix |--> Correlation matrix -- efficiently
#   ## ----------------------------------------------------------------------
#   ## Arguments: V: a covariance matrix (i.e. symmetric and positive definite)
#   ## ----------------------------------------------------------------------
#   ## Author: Martin Maechler, Date: 12 Jun 2003, 11:50
#   p <- (d <- dim(V))[1]
#   if(!is.numeric(V) || length(d) != 2 || p != d[2])
#     stop("`V' is not a square numeric matrix")
#   Is <- sqrt(1/diag(V)) # diag( 1/sigma_i )
#   if(any(!is.finite(Is)))
#     warning("diagonal has non-finite entries")
#   r <- V # keep dimnames
#   r[] <- Is * V * rep(Is, each = p)
#   ##	== D %*% V %*% D  where D = diag(Is)
#   r[cbind(1:p,1:p)] <- 1 # exact in diagonal
#   r
# }
#
# reconst_R <- function(directory,chunksize,chrom){
#   rdsfiles <-dir(directory,pattern = paste0("chr",chrom,".*_",chunksize,".RDS"),full.names = T)
#   maxchunk <-max(c(as.numeric(gsub(paste0(".+chr",chrom,"_([0-9]+)_([0-9]+)_",chunksize,".RDS"),"\\1",rdsfiles))),
#                  as.numeric(gsub(paste0(".+chr",chrom,"_([0-9]+)_([0-9]+)_",chunksize,".RDS"),"\\2",rdsfiles)))
#   diagfiles <- character(maxchunk+1)
#   for(i in 0:maxchunk){
#     diagfiles[i+1] <- rdsfiles[grepl(paste0("chr",chrom,"_",i,"_",i,"_",chunksize,".RDS"),rdsfiles)]
#   }
#   variances<- unlist(lapply(diagfiles,function(x){
#     tdiag <-diag(readRDS(x))
#     return(tdiag[tdiag!=0])
#   }))
#   nSNPs <- nrow(readRDS(rdsfiles[1]))
#   nvars <- list()
#   for(k in 1:length(rdsfiles)){
#     nvars[[as.character(k)]] <- sparse_cov2cor(rdsfiles[k],variances = variances,chunksize = chunksize,nSNPs = nSNPs)
#   }
#   nvars <- lapply(rdsfiles,sparse_cov2cor,variances=variances,chunksize=chunksize,nSNPs=nSNPs)
#   newR <- Reduce("+",nvars)
# }



.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(print(object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

gigsize <- function(rows,cols,size=8){
  (rows/1e9*cols*size)
}



#
# mapfilt <- function(mapfile){
#   mapdf <- read.table(mapfile,header=F,stringsAsFactors = F)
#   mapdf <- rename(mapdf,rsid=V1,pos=V2,map=V3)
#   mapdf <- filter(mapdf,substr(rsid,1,2)=="rs")
#   mapdf <- mutate(mapdf,nrsid=as.numeric(substring(rsid,3)))
# }





