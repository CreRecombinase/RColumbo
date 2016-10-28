#SNPRelate functions

write_gtex_gdsn <- function(geno_txt,geno_gdsn,chunksize=1e5){
  require(readr)
  require(gdsfmt)
  require(tidyr)
  require(dplyr)
  wbrows <- as.integer(strsplit(system(paste0("wc -l ",geno_txt),intern = T),split = " ")[[1]][1])-1
  tgenodat <-   geno_dat <- read_delim(geno_txt,col_names=T,delim="\t",n_max=10)
  sample.id <- colnames(tgenodat)[-1]

  newfile <- createfn.gds(geno_gdsn)
  put.attr.gdsn(newfile$root,"FileFormat","SNP_ARRAY")
  add.gdsn(newfile,"sample.id",sample.id)
  var.id <-add.gdsn(newfile,"snp.id",storage = "uint64",compress = "LZ4_RA:256K")
  var.pos <- add.gdsn(newfile,"snp.position",storage = "uint64",compress="LZ4_RA:256K")
  var.chrom <-add.gdsn(newfile,"snp.chromosome",storage = "uint8",compress="LZ4_RA:256K")
  var.allele <-add.gdsn(newfile,"snp.allele",storage = "string",compress="LZ4_RA:256K")
  var.geno <- add.gdsn(newfile,"genotype",valdim=c(length(sample.id),0),storage="float",compress="LZ4_RA:256K")
  put.attr.gdsn(var.geno,"sample.order")

  gtex_cbf <- function(x,pos){
    annocols <-c("chrom","pos","ref","alt","b37")
    nx <-separate(x,Id,into = annocols,sep="_")
    annodf <- select(nx,one_of(annocols)) %>% mutate(allele=paste0(ref,"/",alt)) %>% select(-ref,-alt,-b37) %>%
      mutate(chrom=ifelse(chrom=="X",23,chrom)) %>% mutate(chrom=as.integer(chrom),pos=as.integer(pos))
    rmna <- is.na(annodf$chrom)|is.na(annodf$pos)
    datamat <- data.matrix(select(nx,-one_of(annocols)))
    datamat <- t(datamat[!rmna,])
    annodf <- filter(annodf,!rmna)
    append.gdsn(var.pos,annodf$pos)
    append.gdsn(var.chrom,annodf$chrom)
    append.gdsn(var.allele,annodf$allele)
    append.gdsn(var.id,(0:(length(annodf$allele)-1))+pos)
    append.gdsn(var.geno,datamat)
  }
  wbdf <- read_delim_chunked(geno_txt,delim = "\t",
                             col_names=T,
                             callback = SideEffectChunkCallback$new(gtex_cbf),
                             chunk_size = chunksize)
  readmode.gdsn(var.pos)
  readmode.gdsn(var.chrom)
  readmode.gdsn(var.allele)
  readmode.gdsn(var.id)
  readmode.gdsn(var.geno)
  closefn.gds(newfile)
  return(T)
}
