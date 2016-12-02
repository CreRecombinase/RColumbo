read_gtex_snp <- function(ngzfile,h5file,chunksize=100000,FlipAllele=T){
  require(readr)
  require(tidyr)
  require(dplyr)
  require(data.table)
  cat("Counting number of rows\n")

  all_genodata <- fread(input = paste0("zcat ",ngzfile),sep = "\t",header = T,data.table = F)
  flipA=FlipAllele
  annocols <-c("chrom","pos","ref","alt","b37")
  all_genodata <-separate(all_genodata,Id,into = annocols,sep="_") %>% mutate(doFlip=as.integer(ref<alt))
  annodf <- select(all_genodata,one_of(annocols),doFlip) %>% select(-ref,-alt,-b37) %>%
    mutate(chrom=ifelse(chrom=="X",23,chrom)) %>% mutate(chrom=as.integer(chrom),pos=as.integer(pos))

  write_h5_df(df = annodf,group = "SNPinfo",outfile = h5file)
  rm(annodf)
  gc()
  all_genodata <-data.matrix(select(all_genodata,-one_of(annocols),-doFlip))
  gc()
  write_dmatrix_h5(h5file = h5file,groupname="SNPdata",dataname = "genotype",Nsnps = nrow(all_genodata),Nind = ncol(all_genodata),data = all_genodata,deflate_level = 4)
  rm(all_genodata)
  return(T)
}


read_gtex_expression <- function(gzfile,h5file){
  require(readr)
  wbrows <- as.integer(strsplit(system(paste0("wc -l ",gzfile),intern = T),split = " ")[[1]][1])
  wbdf <- read_delim(gzfile,delim="\t",col_names = T,guess_max=wbrows)
  colnames(wbdf)[1] <- c("chrom")
  annocols <- c("chrom","start","end","gene_id")
  wbanno <- select(wbdf,one_of(annocols)) %>% mutate(fgeneid=as.integer(gsub("ENSG([0-9]+)\\.[0-9]+","\\1",gene_id))) %>% select(-gene_id) %>%
    mutate(chrom=ifelse(chrom=="X",23,chrom)) %>% mutate(chrom=as.integer(chrom))
  wbexp <- data.matrix(select(wbdf,-one_of(annocols)))
  write_h5_df(wbanno,group = "EXPinfo",outfile = h5file,deflate_level = 4)
  write_dmatrix_h5(h5file = h5file,groupname = "EXPdata",dataname = "expression",Nsnps = nrow(wbexp),Nind = ncol(wbexp),data = wbexp,deflate_level = 4)
  gc()
  return(T)
}
