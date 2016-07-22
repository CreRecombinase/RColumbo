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





