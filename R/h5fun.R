chunk_h5file <- function(in_h5file,out_h5dir,chunksize=NULL){
  require(h5)
  stopifnot(file.exists(in_h5file))

}

exists_group <- function(h5filename,groupname){
  require(h5)
  hff <- h5file(h5filename,'r')
  exg <- existsGroup(hff,groupname = groupname)
  h5close(hff)
  return(exg)
}


set_attr <- function(h5filename,groupname,attrname,attrvalue){
  require(h5)
  hff <- h5file(h5filename,'a')
  h5g <- hff[groupname]
  h5::createAttribute(.Object = h5g,attributename =attrname,data=attrvalue)
  h5close(hff)
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
write_covar_h5 <- function(covarf,h5filename,chunksize=1,deflate_level=9){
  library(dplyr)
  covardat <- read.table(covarf,header=T,stringsAsFactors = F)
  matdat <-t(data.matrix(select(covardat,-ID)))
  write_2dmat_h5(h5filename,"Covardat","covariates",data = matdat)
  h5f <- h5file(h5filename,'a')
  h5g <- createGroup(h5f,"Covarinfo")
  h5g['id'] <- covardat$ID
  # h5createGroup(file = h5file,group = "Covarinfo")
  # h5createDataset(h5file,"/Covarinfo/id",
  #                 dims=c(nrow(covardat)),
  #                 storage.mode="character",size=max(nchar(covardat$ID))+1,
  #                 chunk=chunksize,level=deflate_level)
  # h5write(covardat$ID,file=h5file,name="/Covarinfo/id")
  h5close(h5f)

}


read_ind_h5 <- function(h5filename,groupname,dataname,index){
  require(h5)
  minind <- min(index)
  maxind <- max(index)
  hf <- h5file(h5filename,'r')
  hd <- hf[paste0('/',groupname,'/',dataname)]
  dimd <- dim(hd)
  subd <- hd[,minind:maxind][,(index-minind)+1,drop=F]
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


read_2d_mat_h5 <- function(h5f,groupn,datan,rows=NULL,cols=NULL){
  require(h5)
  stopifnot(file.exists(h5f))
  h5ff <- h5file(h5f,'r')
  h5g <- h5ff[groupn]

  if(is.null(rows)&is.null(cols)){
    data <- h5g[datan][,]
  }else{
    if(is.null(cols)){
      data <- h5g[datan][rows,]
    }else{
      if(is.null(rows)){
        data <- h5g[datan][,cols]
      }else{
        data <- h5g[datan][rows,cols]
      }
    }
  }
  return(data)
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

# write_rows_gzfile <- function(gzfilename,subset_rows,h5filename,groupname,dataname,max.size=335544320L){
#   library(data.table)
#   tempf <- tempfile()
#   write.table(subset_rows,file=tempf,sep="\n",row.names=F,col.names=F,quote=F)
#   cmdtxt <-paste0('awk \'FNR==NR{line[$0]=$0; next} FNR in line\' ',tempf, ' <(gzip -cd ', gzfilename,')')
#   tcmd <- system(cmdtxt,intern = T)
#   all_data <- fread(input = cmdtxt)
#
#   tsh <- system()
#
# }

write_rows_gzfile_h5 <- function(gzfilename,subset_rows,h5filename,groupname,dataname,max.size=335544320L){
  library(iotools)
  gzf <- gzfile(gzfilename,open = 'rt')
  chunk_mat <- iotools::chunk.reader(gzf,sep=" ")
  i <- 1
  mchunk <- iotools::mstrsplit(iotools::read.chunk(chunk_mat,max.size = max.size),sep=" ",type="numeric")
  while(nrow(mchunk)>0){
    cat("chunk:",i,"\n")
    file_rows <- i:(i+nrow(mchunk)-1)
    stopifnot(length(file_rows)==nrow(mchunk))
    mat_rows <- 1:nrow(mchunk)
    sub_chunk <- mat_rows[which(file_rows %in% subset_rows)]
    tchunk <- t(mchunk[sub_chunk,,drop=F])
    cat("dim:",dim(tchunk))
    if(ncol(tchunk)>0){
      write_2dmat_h5(h5filename,groupn=groupname,datan=dataname,data=tchunk,append=T)
    }
    if(i==1){
      nid <- nrow(tchunk)
    }
    i <- file_rows[length(file_rows)]+1
    mchunk <- iotools::mstrsplit(iotools::read.chunk(chunk_mat,max.size = max.size),sep=" ",type="numeric")
  }
  h5f <- h5file(h5filename,'r')
  h5d <- h5f[paste0("/",groupname,"/",dataname)]
  stopifnot(ncol(h5d)==length(subset_rows))
  h5close(h5f)
  return(T)
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






read_df_h5 <- function(h5filepath,groupname=NULL,subcols=NULL,filtervec=NULL){
  require(h5)
  require(dplyr)
  require(lazyeval)
  stopifnot(file.exists(h5filepath))

  f <- h5file(h5filepath,mode = 'r')
  if(!is.null(groupname)){
    stopifnot(existsGroup(f,groupname))
    group <- openGroup(f,groupname)
    dsets <- list.datasets(group)
  }else{
    dsets <- list.datasets(f)
  }
  if(!is.null(subcols)){
    if(!is.null(groupname)){
      subcols <- paste0("/",groupname,"/",subcols)
    }else{
      subcols <- paste0("/",subcols)
    }
    dsets <- dsets[dsets %in% subcols]
  }
  stopifnot(length(dsets)>0)
  if(!is.null(groupname)){
    dsnames <- gsub(pattern = paste0("/",groupname,"/"),"",dsets)
  }else{
    dsnames <-gsub("/","",dsets)
  }
  return(as_data_frame(setNames(lapply(dsets,function(x,file,fvec){
    if(is.null(fvec)){
      return(x=c(file[x][]))
    }else{
      if(fvec[length(fvec)]==(fvec[1]+length(fvec)-1)){
      return(file[x][fvec])
      }else{
        return(file[x][][fvec])
      }
    }
  },file=f,fvec=filtervec),dsnames)
  ))
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




#' read a sparse matrix from HDF5
#' read a sparse matrix from HDF5 stored in Compressed Column Storage (CCS) format
#' @template h5fun
#' @export
read_ccs_h5 <- function(h5filename,groupname,dataname="data",iname="ir",pname="jc"){
  require("h5")
  require("Matrix")
  h5f <- h5::h5file(h5filename,mode = 'r')
  h5g <- h5f[groupname]
  h5attr(h5g,"Layout")
  isSymmetric <- h5attr(h5g,"isSymmetric")=="TRUE"
  data <- h5g[dataname][]
  i <- h5g[iname][]
  p <- h5g[pname][]
  h5::h5close(h5f)
  return(Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,symmetric = isSymmetric))
}

read_coo_h5 <- function(h5filename,groupname,dataname="data",iname="i",jname="j",giveCsparse=F){
  require(h5)
  require(Matrix)
  h5f <- h5::h5file(h5filename,mode = 'r')
  h5g <- h5f[paste0("/",groupname)]
  data_dim <- h5attr(h5g,'Dimensions')
  data <- h5g[dataname][]
  i <- h5g[iname][]
  j <- h5g[jname][]
  h5::h5close(h5f)
  return(Matrix::sparseMatrix(i=i+1,j=j+1,x=data,dims = data_dim,giveCsparse = giveCsparse))
}


append_h5_h5 <- function(input_h5,input_groupname,input_dataname,output_h5,output_groupname,output_dataname){
  library(h5)
  require(Matrix)
  stopifnot(file.exists(input_h5))
  cat("Opening input file\n",input_h5,"\n")
  ih5f <- h5file(name = input_h5,mode = 'r')
  cat("Opening input group\n")
  ih5g <- ih5f[input_groupname]

  ih5d <- ih5g[input_dataname]

  oh5f <- h5file(output_h5,mode='a')
  if(!existsGroup(oh5f,output_groupname)){
    oh5g <- createGroup(oh5f,output_groupname)
  }
  oh5g <- oh5f[output_groupname]
  if(!existsDataSet(.Object = oh5g,datasetname = output_dataname)){
    h5::createDataSet(.Object = oh5g,
                            datasetname = output_dataname,
                            data=ih5d[],
                            chunksize=as.integer(dim(ih5d)/10),
                            maxdimensions = NA_integer_,
                            compression=as.integer(4))
  }else{
    oh5d <- oh5g[output_dataname]
    noh5d <- c(oh5d,ih5d[])
  }
  h5close(ih5f)
  h5close(oh5f)
  gc()
  return(T)
}

write_coo_h5 <- function(h5filename,spmat,groupname,dataname="data",iname="i",jname="j",compression_level=4){
  require(h5)
  spmat <- as(spmat,"dgTMatrix")
  i <- spmat@i
  j <- spmat@j
  data <- spmat@x
  h5f <- h5::h5file(h5filename,'a')
  if(existsGroup(h5f,groupname)){
    h5g <- h5f[groupname]

    tid <- h5g[iname]
    oli <- dim(tid)
    ntidata <- extendDataSet(tid,oli+length(i))
    ntidata[(oli+1):(oli+length(i))] <- i

    tjd <- h5g[jname]
    olj <- dim(tjd)
    ntjdata <- extendDataSet(tjd,olj+length(j))
    ntjdata[(olj+1):(olj+length(j))] <- j

    txd <- h5g[dataname]
    olx <- dim(txd)
    ntxdata <- extendDataSet(txd,olx+length(data))
    ntxdata[(olx+1):(olx+length(data))] <- data
  }else{
    h5g <- createGroup(h5f,groupname)
    #Matrix::sparseMatrix(i=i+1,p=p,x=data))

    h5::createAttribute(.Object = h5g,attributename = "Layout",data="COO")
    h5::createAttribute(.Object = h5g,attributename = "Dimensions",data=dim(spmat))
    id <- h5::createDataSet(.Object = h5g,
                            datasetname = iname,
                            data=i,
                            chunksize=as.integer(length(i)/10),
                            maxdimensions = NA_integer_,
                            compression=as.integer(compression_level))
    pd <- h5::createDataSet(.Object = h5g,
                            datasetname = jname,
                            data=j,
                            chunksize=as.integer(length(j)/10),
                            maxdimensions = NA_integer_,
                            compression=as.integer(compression_level))

    dd <- h5::createDataSet(.Object = h5g,
                            datasetname = dataname,
                            data=data,
                            chunksize=as.integer(length(data)/10),
                            maxdimensions = NA_integer_,
                            compression=as.integer(compression_level))
  }
  h5close(h5f)
  return(T)


}

write_ccs_h5 <- function(h5filename,spmat,groupname,dataname="data",iname="ir",pname="jc",compression_level=8,symmetric=F){
  require(h5)
  h5f <- h5::h5file(h5filename,'a')
  h5g <- createGroup(h5f,groupname)
  #Matrix::sparseMatrix(i=i+1,p=p,x=data))
  i <- spmat@i
  p <- spmat@p
  data <- spmat@x
  h5::createAttribute(.Object = h5g,attributename = "Layout",data="CCS")
  h5::createAttribute(.Object = h5g,attributename = "Dimensions",data=dim(spmat))
  h5::createAttribute(.Object = h5g,attributename = "isSymmetric",data=ifelse(symmetric,"TRUE","FALSE"))
  id <- h5::createDataSet(.Object = h5g,
                          datasetname = iname,
                          data=i,
                          chunksize=as.integer(length(i)/10),
                          maxdimensions = NA_integer_,
                          compression=as.integer(compression_level))
  pd <- h5::createDataSet(.Object = h5g,
                          datasetname = pname,
                          data=p,
                          chunksize=as.integer(length(p)/10),
                          maxdimensions = NA_integer_,
                          compression=as.integer(compression_level))

  dd <- h5::createDataSet(.Object = h5g,
                          datasetname = dataname,
                          data=data,
                          chunksize=as.integer(length(data)/10),
                          maxdimensions = NA_integer_,
                          compression=as.integer(compression_level))

  h5close(h5f)
  return(T)
}


read_attr <- function(h5filename,datapath,attrname){
  requireNamespace("h5")
  h5f <- h5::h5file(h5filename,'r')
  datag <- h5f[datapath]
  retval <- h5attr(datag,attrname)
  h5::h5close(h5f)
  return(retval)
}

read_vec <- function(h5filename,datapath){
  requireNamespace("h5")
  h5f <- h5::h5file(h5filename,'r')
  data <- h5f[datapath][]
  h5::h5close(h5f)
  return(data)
}

write_vec <- function(h5filename,vec,datapath,deflate_level=4,chunksize=NULL){
  require(h5)
  h5f <- h5::h5file(h5filename,'a')
  if(is.null(chunksize)){
    chunksize <- length(vec)/2
  }
  data <- h5f[datapath,compression=deflate_level,chunksize=chunksize] <- vec
  h5::h5close(h5f)
}





