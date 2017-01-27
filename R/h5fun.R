chunk_h5file <- function(in_h5file,out_h5dir,chunksize=NULL){
  require(h5)
  stopifnot(file.exists(in_h5file))

}

chunk_h5d <- function(in_h5file,in_groupname,in_dataname,out_h5dir,out_fname,out_groupname,out_dataname,chunksize=NULL){
  require(h5)
  stopifnot(!is.null(chunksize))
  df <- h5file(in_h5file,'r')
  dpath <- df[paste0('/',in_groupname,'/',in_dataname)]
  ddim <- dim(dpath)

}




sample_h5_df <-function(mh5file,groupname,blocksize=20000,pval.cutoff=1e-5,gwasdf=NULL,tstat_df=337){
  require(rhdf5)
  cut_t <- qt(p = pval.cutoff,df = tstat_df,lower.tail = F)
  tdf <- read_h5_df(mh5file,groupname=groupname,subcols = c("theta","serr"))

  filtervec <- abs(tdf$theta/tdf$serr)>cut_t
  rm(tdf)
  snpdf <- read_h5_df(mh5file,groupname=groupname,subcols=c("chrom","pos")) %>% distinct
  stopifnot(nrow(distinct(snpdf,chrom))==1)
  ldblocks <- sort(unique(c(seq(from=min(snpdf$pos)-1,to=max(snpdf$pos)+1,by = blocksize),max(snpdf$pos)+1)))
  snpdf <- mutate(snpdf,ldblock=cut(pos,breaks = ldblocks,labels = F,right=F))
  sigdf <- read_h5_df(mh5file,groupname=groupname,filtervec = filtervec)  %>% inner_join(snpdf,by=c("chrom","pos"))
  if(!is.null(gwasdf)){
    sigdf <- inner_join(sigdf,gwasdf,by=c("chrom","pos"))
  }
  sigdf<-  mutate(sigdf,tstat=theta/serr) %>% mutate(pvale=pt(abs(tstat),tstat_df,lower.tail = F)) %>%
    group_by(fgeneid,ldblock) %>% filter(pvale==min(pvale)) %>% ungroup()
  gc()
  return(sigdf)
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


read_ind_h5 <- function(h5filename,groupname,dataname,index){
  require(h5)
  minind <- min(index)
  maxind <- max(index)
  hf <- h5file(h5filename,'r')
  hd <- hf[paste0('/',groupname,'/',dataname)]
  dimd <- dim(hd)
  subd <- hd[,minind:maxind][,(index-minind)+1]
  h5close(hf)
  return(subd)
}



read_subset_h5 <- function(h5filename,groupname,dataname,index,chunk_size=NULL,chunk_num=NULL){
  stopifnot(!is.null(chunk_size)||!is.null(chunk_num))
  require(BBmisc)
  if(!is.null(chunk_size)){
    chunkind <-chunk(x = index,chunk.size = chunk_size)
  }else{
    chunkind <-chunk(x = index,n.chunks = chunk_num)
  }
  fdata <- do.call("cbind",
                   lapply(chunkind,read_ind_h5,h5filename=h5filename,groupname=groupname,dataname=dataname))
  stopifnot(length(index)==ncol(fdata))
  return(fdata)
}

copy_subset_h5 <-function(h5filename,groupname,dataname,index,chunk_size=NULL,chunk_num=NULL,newh5filename){
  stopifnot(!is.null(chunk_size)||!is.null(chunk_num))
  require(BBmisc)
  if(!is.null(chunk_size)){
    chunkind <-chunk(x = index,chunk.size = chunk_size)
  }else{
    chunkind <-chunk(x = index,n.chunks = chunk_num)
  }
  for(i in 1:length(chunkind)){
    tdata <- read_ind_h5(h5filename,groupname,dataname,chunkind[[i]])
    write_2dmat_h5(newh5filename,groupname,dataname,deflate_level=4,data=tdata,append=T)
  }
  nf <- h5file(newh5filename,mode='r')
  nd <- nf[paste0('/',groupname,'/',dataname)]
  stopifnot(ncol(nd)==length(index))
  h5close(nf)
  return(T)
}



write_2dmat_h5 <- function(h5f,groupn,datan,chunksize=c(5000,5000),deflate_level=4,data=NULL,append=F){
  require(h5)
  stopifnot(!is.null(data))
  if(any(chunksize>dim(data))){
    chunksize[which(chunksize>=dim(data))] <-dim(data)[which(chunksize>=dim(data))]/2
  }
  hf <- h5file(h5f,mode='a')
  if(!existsGroup(hf,groupname=groupn)){
    grp <- createGroup(hf,groupname = groupn)
  }else{
    grp <- openGroup(hf,groupname = paste0("/",groupn))
  }
  if(existsDataSet(grp,datan)){
    if(append){
      ds <- openDataSet(grp,datan)
      current_row <- nrow(ds)
      current_col <- ncol(ds)

      data_row <- nrow(data)
      data_col <- ncol(data)
      stopifnot(current_row==data_row)
      new_dim <- c(current_row,current_col+data_col)
      nds <-extendDataSet(ds,dims = new_dim)
      nds[,(current_col+1):new_dim[2]] <- data
    }
    else{stop("Dataset already exists!")}
  }else{
  wdata <- createDataSet(grp,datasetname = datan,data=data,chunksize=chunksize,compression=deflate_level)
stopifnot(all(dim(wdata)==dim(data)))
  }
  h5close(hf)
  return(T)
}


write_covar_2d <- function(covarf,h5filename,chunksize=1,deflate_level=9){
  require(h5)
  covardat <- read.table(covarf,header=T,stringsAsFactors = F)
  matdat <-t(data.matrix(select(covardat,-ID)))
  ncovar <- ncol(matdat)
  nid <- nrow(matdat)
  write_2dmat_h5(h5f = h5filename,groupn = "Covardat",datan = "covariates",data=matdat,append=F)
  h5createGroup(file = h5file,group = "Covarinfo")
  h5createDataset(h5file,"/Covarinfo/id",
                  dims=c(nrow(covardat)),
                  storage.mode="character",size=max(nchar(covardat$ID))+1,
                  chunk=chunksize,level=deflate_level)
  h5write(covardat$ID,file=h5file,name="/Covarinfo/id")
}


 write_gtex_eqtl_h5 <- function(txtfile,eqtlh5file,chunksize=500000){
   require(readr)
   require(tidyr)
   require(dplyr)


   gtex_cbf <- function(x,pos){
     th5file <- eqtlh5file
     annocols <-c("chrom","pos","ref","alt","b37")
     nx <-separate(x,variant_id,into = annocols,sep="_",convert = T) %>% mutate(doFlip=as.integer(ref<alt),fgeneid=as.integer(gsub("ENSG([0-9]+).+","\\1",gene_id)))
     annodf <- select(nx,-b37,-gene_id,-ref,-alt) %>%rename(pval.e=pval_nominal,weight=slope) %>%
       mutate(chrom=ifelse(chrom=="X",23,chrom)) %>%
       mutate(chrom=as.integer(chrom),
              pos=as.integer(pos),
              weight=ifelse(doFlip==1,-weight,weight))
     rmna <- is.na(annodf$chrom)|is.na(annodf$pos)
     annodf <- filter(annodf,!rmna)
     write_df_h5(df = annodf,groupname = "cis_eQTL",outfile = th5file,deflate_level = 4)
   }
   wbdf <- read_delim_chunked(txtfile,delim = "\t",col_names=T,callback = SideEffectChunkCallback$new(gtex_cbf),chunk_size = chunksize)
   gc()
   return(T)
 }


 write_df_h5 <- function(df,groupname,outfile,deflate_level=4L,chunksize=1000L){
   if(nrow(df)<chunksize){
     chunksize <- nrow(df)
   }
   deflate_level <- as.integer(deflate_level)
   require(h5)
   dataname <- colnames(df)
   f <-h5file(outfile,mode = 'a')
   if(existsGroup(f,groupname)){
     h5close(f)
     res <- append_df_h5(df,groupname,outfile,deflate_level)
     return(res)
   }
   group <- createGroup(f,groupname)
   for(i in 1:length(dataname)){
     dsn <- dataname[i]
     td <- df[[dsn]]
     tdata <- createDataSet(.Object = group,
                            datasetname = dsn,
                            type = typeof(td),
                            dimensions = length(td),
                            chunksize =chunksize,
                            maxdimensions = NA_integer_,
                            compression = deflate_level)
     tdata[] <- td
   }
   h5close(f)
   return(T)
 }

append_df_h5 <- function(df,groupname,outfile,deflate_level=4){
  require(h5)
  dataname <- colnames(df)
  f <-h5file(outfile,mode = 'a')
  group <- openGroup(f,groupname)
  for(i in 1:length(dataname)){
    dsn <- dataname[i]
    td <- df[[dsn]]
    if(existsDataSet(group,dsn)){
      tdata <- openDataSet(group,dsn)
      odim <- dim(tdata)
      ntdata <- extendDataSet(tdata,odim+length(td))
      ntdata[(odim+1):(odim+length(td))] <- td
    }else{
      h5close(f)
      stop(paste0(dsn," does not exist in ",outfile))
    }
  }
  h5close(f)
  return(T)
}






read_df_h5 <- function(h5filepath,groupname,subcols=NULL,filtervec=NULL){
  require(h5)
  require(dplyr)
  require(lazyeval)
  stopifnot(file.exists(h5filepath))

  f <- h5file(h5filepath,mode = 'r')
  stopifnot(existsGroup(f,groupname))
  group <- openGroup(f,groupname)
  dsets <- list.datasets(group)
  if(!is.null(subcols)){
    subcols <- paste0("/",groupname,"/",subcols)
    dsets <- dsets[dsets %in% subcols]
  }
  stopifnot(length(dsets)>0)
  dsnames <- gsub(pattern = paste0("/",groupname,"/"),"",dsets)
  return(as_data_frame(setNames(lapply(dsets,function(x,file,fvec){
    if(is.null(fvec)){
      return(x=file[x][])
    }else{
      return(file[x][filtervec])
    }
  },file=f,fvec=filtervec),dsnames)
  ))
}


read_h5_df <- function(h5file,groupname,subcols=NULL,filtervec=NULL){
  require(rhdf5)
  require(dplyr)
  require(purrr)
  stopifnot(file.exists(h5file))
  h5file <- system(paste0("echo ",h5file),intern = T)
  datarows <- h5ls(h5file) %>% filter(group==paste0("/",groupname)) %>% mutate(h5file=h5file,dim=as.integer(dim)) %>%
    filter(dclass %in% c("INTEGER","FLOAT","DOUBLE"))
  if(!is.null(subcols)){
    datarows <- filter(datarows,name %in% subcols)
  }
  retdf <- datarows %>% split(.$name) %>% map(~read_h5_type(h5file=h5file,groupname=groupname,name=.$name,dclass=.$dclass,filtervec=filtervec)) %>% bind_cols()
  return(retdf)
}







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





