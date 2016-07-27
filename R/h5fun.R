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


#' Write haplotype data to HDF5 data
#' @param haplotype matrix (one SNP per column, one individual per row)
#' @param legdat legend dataframe
#' @param hap_h5file output HDF5 filename
#' @param chunksize chunk size for compression
write_hap_h5 <- function(haplotype,legdat,hap_h5file,chunksize=1000){
  if(nrow(haplotype)==nrow(legdat)){
    haplotype <- t(haplotype)
  }
  chunkvec <-c(nrow(haplotype),chunksize)
  h5createFile(hap_h5file)
  maxchar <- max(nchar(legdat$ID))
  h5createGroup(file = hap_h5file,group = "Genotype")
  h5createGroup(file = hap_h5file,group = "Legend")
  h5createDataset(hap_h5file,"/Genotype/Haplotype",
                  dims = c(nrow(haplotype),ncol(haplotype)),
                  storage.mode = "integer",
                  chunk=chunkvec,level=2)
  h5createDataset(hap_h5file,"/Legend/rsid",
                  dims=c(nrow(legdat)),
                  storage.mode="character",
                  size=max(nchar(legdat$ID)),
                  chunk=chunksize,level=2)
  h5createDataset(hap_h5file,"/Legend/pos",
                  dims=c(nrow(legdat)),
                  storage.mode="integer",
                  size=maxchar,
                  chunk=chunksize,level=2)
  h5createDataset(hap_h5file,"/Legend/allele0",
                  dims=c(nrow(legdat)),
                  storage.mode="character",
                  size=max(nchar(legdat$allele0)),
                  chunk=chunksize,level=2)
  h5createDataset(hap_h5file,"/Legend/allele1",
                  dims=c(nrow(legdat)),
                  storage.mode="character",
                  size=max(nchar(legdat$allele1)),
                  chunk=chunksize,level=2)

  h5write(legdat$ID,file=hap_h5file,name="/Legend/rsid")
  h5write(legdat$pos,file=hap_h5file,name="/Legend/pos")
  h5write(legdat$allele0,file=hap_h5file,name="/Legend/allele0")
  h5write(legdat$allele1,file=hap_h5file,name="/Legend/allele1")
  h5write(haplotype,file=hap_h5file,name="/Genotype/Haplotype")
}

write_file_h5 <- function(hapfile,legfile,hap_h5file,chunksize=1000){
  require(rhdf5)
  hapdata <- data.matrix(read_delim(hapfile,col_names=F,delim=" "))
  legdata <- read_delim(legendfile,col_names=T,delim=" ")
  hapdata <- hapdata[!duplicated(legdata$ID),]
  legdata <- legdata[!duplicated(legdata$ID),]
  write_hap_h5(haplotype = hapdata,hap_h5file = hap_h5file,legdat = legdata,chunksize = chunksize)
}


read_h5_df <- function(h5file,group){
  td <-h5read(h5file,group)
  retdf <- bind_cols(mapply(function(x,y){
    ret <- data_frame(c(x))
    colnames(ret)=y
    return(ret)
  },td,names(td),USE.NAMES=F))
  return(retdf)
}


#' Write haplotype data to HDF5 data
#' @param haplotype matrix (one SNP per column, one individual per row)
#' @param h5file file output file
#' @param chunksize
read_hap_h5 <- function(hap_h5file,rslist=NULL,poslist=NULL){
  stopifnot(length(rslist>0)|length(poslist>0))
  legdata <-read_h5_df(hap_h5file,"Legend")
  if(length(rslist)>0){
    subset_ind <- which(legdata$rsid %in% rslist)
    legdata <- filter(legdata,rsid %in% rslist)
  }
  hapdim <-h5ls(hap_h5file) %>% filter(name=="Haplotype")
  nind <- as.numeric(gsub("^([0-9]+) x [0-9]+$","\\1",hapdim$dim))
  shapdat <- h5read(hap_h5file,name="/Genotype/Haplotype")
  shapdat <- t(shapdat[,subset_ind])
  return(shapdat)
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


#
# mapfilt <- function(mapfile){
#   mapdf <- read.table(mapfile,header=F,stringsAsFactors = F)
#   mapdf <- rename(mapdf,rsid=V1,pos=V2,map=V3)
#   mapdf <- filter(mapdf,substr(rsid,1,2)=="rs")
#   mapdf <- mutate(mapdf,nrsid=as.numeric(substring(rsid,3)))
# }





