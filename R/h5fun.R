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

read_h5_df_l <- function(h5files,chroms=1:22,group){
  stopifnot(length(h5files)==length(chroms))
  retdf <- bind_rows(mapply(function(chromi,file,group){
    tdf <- read_h5_df(file,groupname = group) %>%mutate(chrom=chromi,ind=1:n())
    return(tdf)
  },chroms,h5files,MoreArgs=list(group=group),SIMPLIFY = F))
}


write_geno_h5 <- function(mapfile,gengz,famfile,dbsnpfile,h5file,chunksize=100000){
  require(dplyr)
  legdf <- read.table(mapfile,header=F,stringsAsFactors = F) %>% rename(chrom=V1,rsid=V2,score=V3,pos=V4) %>% mutate(rsid=as.integer(substring(rsid,3)))
  legdf <- select(legdf,-score)
  nsnp <- nrow(legdf)
  famdf <- read.table(famfile,header=T,stringsAsFactors=F) %>% slice(-1)
  nind <- nrow(famdf)
  write_genofile_h5(gengz,dbsnpfile,nind,nsnp,chunksize,h5file,2)
  write_h5_df(df = legdf,group = "Legend",outfile = h5file,deflate_level = 3)
  nlegdf <- read_h5_df(h5file,"Legend")
  tanno <- read_delim(gengz,delim = " ",col_names = c("chrom","rsid","pos","ref","alt"))
}






sample_h5_df <-function(mh5file,groupname,blocksize=20000,pval.cutoff=1e-5,gwasdf=NULL){
  require(rhdf5)
  cut_t <- qt(p = pval.cutoff,df = 337,lower.tail = F)
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
  sigdf<-  mutate(sigdf,tstat=theta/serr) %>% mutate(pvale=pt(abs(tstat),337,lower.tail = F)) %>%
    group_by(fgeneid,ldblock) %>% filter(pvale==min(pvale)) %>% ungroup()
  gc()
  return(sigdf)
}



fsample_h5_df <-function(mh5file,groupname,blocksize=20000,pval.cutoff=1e-5,gwasdf=NULL,chunksize=500000){
  require(rhdf5)
  cut_t <- qt(p = pval.cutoff,df = 337,lower.tail = F)
  tdf <- fread_h5_df(mh5file,groupname=groupname,subcols = c("theta","serr"),chunksize=chunksize)
  filtervec <- abs(tdf$theta/tdf$serr)>cut_t
  filvi <- which(filtervec)

  rm(tdf)
  snpdf <- fread_h5_df(mh5file,groupname=groupname,subcols=c("chrom","pos"),chunksize=chunksize) %>% distinct
  stopifnot(nrow(distinct(snpdf,chrom))==1)
  cat(snpdf$chrom[1],"\n")
  ldblocks <- sort(unique(c(seq(from=min(snpdf$pos)-1,to=max(snpdf$pos)+1,by = blocksize),max(snpdf$pos)+1)))
  snpdf <- mutate(snpdf,ldblock=cut(pos,breaks = ldblocks,labels = F,right=F))
  # sigdf <- read_h5_df(mh5file,groupname=groupname,filtervec = filtervec) %>% inner_join(snpdf,by=c("chrom","pos"))
  sigdf <- fread_h5_df(mh5file,groupname=groupname,indexvec = filvi,chunksize=chunksize)%>% inner_join(snpdf,by=c("chrom","pos"))
  if(!is.null(gwasdf)){
    sigdf <- inner_join(sigdf,gwasdf,by=c("chrom","pos"))
  }
  sigdf<-  mutate(sigdf,tstat=theta/serr) %>% mutate(pvale=pt(abs(tstat),337,lower.tail = F)) %>%
    group_by(fgeneid,ldblock) %>% filter(pvale==min(pvale)) %>% ungroup()
  gc()
  return(sigdf)
}



# construct_LD_graph <- function(geno_h5,geno_leg,LD_cutoff,chunksize=30000){
#   require(graph)
#   require(RBGL)
#   require(tidyr)
#   require(SparseM)
#
#
#   stopifnot(all(c("chrom","pos","snp_ind") %in% colnames(geno_leg)))
#   stopifnot(all(!duplicated(geno_leg$snp_ind)))
#   geno_leg <- group_by(geno_leg,chrom) %>% mutate(chunkid=cut(snp_ind,breaks = ceiling(n()/chunksize),labels = F)) %>% ungroup()
#
#   dimtot <-nrow(geno_leg)
#   group_by(geno_leg,chrom,chunkid)
#   tgeno_leg <- filter(geno_leg,chrom==1,chunkid==1)
#   tgenod <- read_fmat_chunk_ind(geno_h5,"SNPdata","genotype",tgeno_leg$snp_ind[1:100])
#   genoco <- cor(tgenod)
#   genoco <- as_data_frame(genoco)
#   ogenoco <- as_data_frame(gtex_LDM)
#   colnames(ogenoco) <- paste0("S",1:ncol(ogenoco))
#   ogenoco <- mutate(ogenoco,rowid=1:n())
#   ogenoco <- gather(ogenoco,colid,LD,-rowid)
#   ogenoco <- mutate(ogenoco,colid=as.integer(gsub(pattern = "S([0-9]+)",replacement = "\\1",x = colid)))
#
#
#
#   colnames(genoco) <- paste0("S",1:ncol(genoco))
#   genoco <- mutate(genoco,rowid=1:n())
#   genoco <- gather(genoco,colid,LD,-rowid)
#   genoco <- mutate(genoco,colid=as.integer(gsub(pattern = "S([0-9]+)",replacement = "\\1",x = colid)))
#
#
#   # genoc <- cor_h5(geno_h5,groupname = "SNPdata",dataname = "genotype",indvec = tgeno_leg$snp_ind,LDcutoff = 0.01,cutBelow = T)
#
#
#   tgenoc <- arrange(genoc,rowind+colind) %>%slice(1:10000)
#
# ogenoco%>%ggplot(aes(rowid,colid))+geom_raster(aes(fill=LD))
#   snpa <- select(tgeno_leg,chromA=chrom,posA=pos,rowind=snp_ind)
#   snpb <-select(tgeno_leg,chromB=chrom,posB=pos,colind=snp_ind)
#   ngenoc <- inner_join(genoc,snpa) %>% inner_join(snpb)
#   ngenoc <- filter(ngenoc,chromA==chromB,posA<posB)
#   edgel <- as.matrix.csr(new("matrix.coo",
#                ra=rep.int(1,length(ngenoc$rowind)),
#                ja=ngenoc$colind,
#                ia=ngenoc$rowind,
#                dimension=c(dimtot,dimtot)))
#   genog <- graphAM(adjMat = edgel,
#                    edgemode ="directed")
#
#
#   genog <- sparseM2Graph(edgel,nodeNames =as.character(geno_leg$snp_ind),
#                 edgemode = "directed")
#
#
# }






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



read_geno_h5_pos <- function(h5file,query_posdf){
  #Figure out which dataset has the same dimensions as the array dataset and is named 'pos', use that to index
  require(rhdf5)
  datainfo <- h5ls(h5file)
  arrayrow <- filter(datainfo,dclass=="ARRAY")
  arraygroup <- substring(arrayrow$group,2)
  arraydata <- arrayrow$name
  annorows <- filter(datainfo,dim==arrayrow$dim,name %in% colnames(query_posdf))
  annogroup <- substring(annorows$group[1],2)
  tannodf <- read_h5_df(h5file,annogroup) %>% mutate(annopos=1:n())
  read_ind <- semi_join(tannodf,query_posdf) %>% distinct()
  # read_ind <- read_ind[!duplicated(annopos[read_ind])]
  retmat <- read_dmat_chunk_ind(h5file,arraygroup,arraydata,read_ind$annopos)
  colnames(retmat) <- as.character(read_ind$pos)
  return(retmat)
}





read_gtex_expression <- function(gzfile,h5file){
  require(readr)
  wbrows <- as.integer(strsplit(system(paste0("wc -l ",gzfile),intern = T),split = " ")[[1]][1])
  wbdf <- read_delim(gzfile,delim="\t",col_names = T,guess_max=wbrows)
  colnames(wbdf)[1] <- c("chr")
  annocols <- c("chr","start","end","gene_id")
  wbanno <- select(wbdf,one_of(annocols)) %>% mutate(fgeneid=as.integer(gsub("ENSG([0-9]+)\\.[0-9]+","\\1",gene_id))) %>% select(-gene_id) %>%
    mutate(chr=ifelse(chr=="X",23,chr)) %>% mutate(chr=as.integer(chr))
  wbexp <- data.matrix(select(wbdf,-one_of(annocols)))
  write_h5_df(wbanno,group = "EXPinfo",outfile = h5file,deflate_level = 4)
  write_dmatrix_h5(h5file = h5file,groupname = "EXPdata",dataname = "expression",Nsnps = nrow(wbexp),Nind = ncol(wbexp),data = wbexp,deflate_level = 4)
  gc()
  return(T)
}



read_dgn_exp <- function(ngzfile,h5file,chunksize=100000){




}


read_gtex_snp <- function(ngzfile,h5file,chunksize=100000){
  require(readr)
  require(tidyr)
  require(dplyr)
  wbrows <- as.integer(strsplit(system(paste0("wc -l ",ngzfile),intern = T),split = " ")[[1]][1])
  gtex_cbf <- function(x,pos){
    th5file <- h5file
    annocols <-c("chrom","pos","ref","alt","b37")
    nx <-separate(x,Id,into = annocols,sep="_") %>% mutate(doFlip=as.integer(ref<alt))
    annodf <- select(nx,one_of(annocols),doFlip) %>% select(-ref,-alt,-b37) %>%
      mutate(chrom=ifelse(chrom=="X",23,chrom)) %>% mutate(chrom=as.integer(chrom),pos=as.integer(pos))
    rmna <- is.na(annodf$chrom)|is.na(annodf$pos)
    datamat <- data.matrix(select(nx,-one_of(annocols),-doFlip))
     datamat[nx$doFlip==1,] <-abs(2-datamat[nx$doFlip==1,])
    datamat <- datamat[!rmna,]
    annodf <- filter(annodf,!rmna)
    write_h5_df(df = annodf,group = "SNPinfo",outfile = th5file,deflate_level = 4)
    write_dmatrix_h5(h5file = th5file,groupname="SNPdata",dataname = "genotype",Nsnps = wbrows,Nind = ncol(datamat),data = datamat,deflate_level = 4)
  }
  wbdf <- read_delim_chunked(ngzfile,delim = "\t",col_names=T,callback = SideEffectChunkCallback$new(gtex_cbf),chunk_size = chunksize)
  gc()
  return(T)
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


fread_h5_df <- function(h5file,groupname,subcols=NULL,indexvec=numeric(0),chunksize=500000){
  require(rhdf5)
  require(dplyr)
  require(purrr)
  stopifnot(file.exists(h5file))
  if(length(indexvec)>0){
    if(length(indexvec)>chunksize){
      mbreaks <-c(seq(min(indexvec)-1,max(indexvec),by=chunksize),max(indexvec)+10)
      chunkl <- split(indexvec,f = cut(indexvec,mbreaks,labels = F))
    }else{
      chunkl <- list(indexvec)
    }
  }
  else{
    chunkl <- list(numeric(0))
  }
  h5file <- system(paste0("echo ",h5file),intern = T)
  datarows <- h5ls(h5file) %>% filter(group==paste0("/",groupname)) %>% mutate(h5file=h5file,dim=as.integer(dim)) %>%
    filter(dclass %in% c("INTEGER","FLOAT","DOUBLE"))
  if(!is.null(subcols)){
    datarows <- filter(datarows,name %in% subcols)
  }
  # Rcpp::DataFrame read_h5_df_col(const std::string h5file, const std::string groupname, const std::string dataname,const std::string dclass,const arma::uvec col_index){
  # name <- datarows$name[1]
  # dclass <- datarows$dclass[1]
  # tres <-   read_h5_df_col(h5file=h5file,groupname=groupname,dataname=name,dclass=dclass,col_index=filtervec)
  retdf <- bind_rows(lapply(chunkl,function(x,datar){
    return(datar %>% split(.$name) %>% map(~read_h5_df_col(h5file=h5file,groupname=groupname,dataname=.$name,dclass=.$dclass,col_index=x)) %>% bind_cols())
  },datar=datarows))
  return(retdf)

}


read_h5_type <-function(h5file,groupname,name,dclass,filtervec=NULL){
  require(dplyr)
  # groupname <- substring(group,2)
  if(dclass=="FLOAT"){
    if(!is.null(filtervec)){
      return(read_Rfloat_h5(h5file,groupname,name)[filtervec])
    }else{
      return(read_Rfloat_h5(h5file,groupname,name))
    }
  }else{
    if(dclass=="INTEGER"){
      if(!is.null(filtervec)){
        return(read_Rint_h5(h5file,groupname,name)[filtervec])
      }else{
        return(read_Rint_h5(h5file,groupname,name))
      }
    }else{
      stop(paste0("Type:",dclass," not supported by read_h5_type\n"))
    }
  }
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


read_fgid <- function(res_row){
  library(tidyr)
  snpleg <-read_h5_df(h5file,"SNPinfo") %>% mutate(snp_ind=1:n())
  expleg <-read_h5_df(h5file,"EXPinfo") %>% mutate(exp_ind=1:n()) %>% filter(fgeneid %in% fgeneidv)
  snpexp <- group_by(expleg,fgeneid) %>%
    do(snpdf=filter(snpleg,chrom==.$chr,
                    (abs(pos-.$start)<cisdist|abs(pos-.$end)<cisdist))) %>%
    unnest(snpdf) %>% inner_join(expleg)
  tidydat <- function(df,valname="snp"){
    colnames(df) <- paste0("ind:",1:ncol(df))
    df <- mutate(df,row=paste0(valname,":",1:n()))
    retdf <- gather(df,ind,value,-row)

  }

  df <- as_data_frame(read_fmat_chunk_ind(h5file,"SNPdata","genotype",indvec=snpexp$snp_ind[snpexp$fgeneid==fgeneidv[1]]))
  snpdat <- group_by(snpexp,fgeneid) %>% do(snpdata=as_data_frame(read_fmat_chunk_ind(h5file,"SNPdata","genotype",indvec = .$snp_ind)))
  snpdat <- ungroup(snpdat)
  expdat <- read_fmat_chunk_ind(h5file,"EXPdata","expression",indvec=expleg$exp_ind)
  snpdat <- unnest(snpdat)







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





