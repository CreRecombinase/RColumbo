
read_1kg_gz <- function(data_gzfile,leg_gzfile,map_gzfile,h5filename){
  require(readr)
  require(dplyr)
  require(h5)
  require(iotools)
  chr <- as.integer(gsub(".+EUR.chr([0-9]+)_1kg.+","\\1",data_gzfile))
  leg_data <- read_delim(leg_gzfile,delim=" ",col_names = c("rsid","pos","ref","alt"),skip = 1) %>%
    mutate(snp_id=1:n(),chrom=chr)
  map_data <- read_delim(map_gzfile,delim=" ",col_names = c("rsid","pos","map")) %>% mutate(chrom=chr,map_id=1:n())
  leg_map <- inner_join(leg_data,map_data) %>% distinct(snp_id,.keep_all=T) %>% distinct(map_id,.keep_all=T) %>% distinct(chrom,pos,.keep_all=T)
  w_leg_map <- select(leg_map,rsid,chrom,pos,ref,alt,map)
  write_df_h5(df = w_leg_map,outfile = h5filename,groupname = "SNPinfo",deflate_level = 4)
  gzf <- gzfile(data_gzfile,open = 'rt')
  chunk_mat <- chunk.reader(gzf,sep=" ")
  chunkr <- 0
  i <- 0
  mchunk <- mstrsplit(read.chunk(chunk_mat,max.size = 335544320L),sep=" ",type="numeric")
  while(nrow(mchunk)>0){
    cat("chunk:",i,"\n")
    sub_chunk <- leg_map$snp_id[leg_map$snp_id>chunkr&leg_map$snp_id<=(chunkr+nrow(mchunk))]-chunkr
    tchunk <- t(mchunk[sub_chunk,])
    cat("dim:",dim(tchunk))
    if(ncol(tchunk)>0){
      write_2dmat_h5(h5filename,groupn="SNPdata",datan="genotype",data=tchunk,append=T)
    }
    if(i==0){
      nid <- nrow(tchunk)
    }
    i <- i+1
    chunkr <- chunkr+nrow(mchunk)
    mchunk <- mstrsplit(read.chunk(chunk_mat,max.size = 335544320L),sep=" ",type="numeric")
  }
  return(T)
}

subset_map_gz <- function(map_gzfile,h5filename,newh5filename){
  require(readr)
  require(dplyr)
  chr <- as.integer(gsub(".+chr([0-9]+).interpolated_genetic.+","\\1",map_gzfile))

  snp_leg <- read_df_h5(h5filename,groupname = "SNPinfo")
  new_snp_leg <- inner_join(snp_leg,map_data)
  copy_subset_h5(h5filename = h5filename,groupname = "SNPdata",dataname = "genotype",index = new_snp_leg$snp_id,chunk_num = 5,newh5filename = newh5filename)
  new_snp_leg <- select(new_snp_leg,chrom,pos,ref,alt,map)
}







read_gtex_expression <- function(gzfile,h5filename){
  require(readr)
  require(dplyr)
  require(h5)
  wbrows <- as.integer(strsplit(system(paste0("wc -l ",gzfile),intern = T),split = " ")[[1]][1])
  wbdf <- read_delim(gzfile,delim="\t",col_names = T)
  colnames(wbdf)[1] <- c("chrom")
  annocols <- c("chrom","start","end","gene_id")
  wbanno <- select(wbdf,one_of(annocols)) %>%
    mutate(fgeneid=as.integer(gsub("ENSG([0-9]+)\\.[0-9]+","\\1",gene_id))) %>%
    select(-gene_id) %>%
    mutate(chrom=ifelse(chrom=="X",23,chrom)) %>%
    mutate(chrom=as.integer(chrom))
  wbexp <- t(data.matrix(select(wbdf,-one_of(annocols))))
  write_df_h5(wbanno,group = "EXPinfo",outfile = h5filename,deflate_level = 4)

  write_dmatrix_h5(h5file = h5filename,groupname = "EXPdata",dataname = "expression",Nsnps = ncol(wbexp),Nind = nrow(wbexp),data = wbexp,deflate_level = 4)

  sexp <- read_dmat_chunk_ind(h5filename,groupname="EXPdata",dataname="expression",indvec = 1:10)

  return(TRUE)
}


read_gtex_expression_2d <- function(gzfile,h5filename){
  require(readr)
  require(dplyr)
  require(h5)
  gzfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Analysis_v6p_eQTL_expression_matrices/Muscle_Skeletal_Analysis.v6p.normalized.expression.bed.gz"
  wbrows <- as.integer(strsplit(system(paste0("wc -l ",gzfile),intern = T),split = " ")[[1]][1])
  wbdf <- read_delim(gzfile,delim="\t",col_names = T)
  colnames(wbdf)[1] <- c("chrom")
  annocols <- c("chrom","start","end","gene_id")
  wbanno <- select(wbdf,one_of(annocols)) %>%
    mutate(fgeneid=as.integer(gsub("ENSG([0-9]+)\\.[0-9]+","\\1",gene_id))) %>%
    select(-gene_id) %>%
    mutate(chrom=ifelse(chrom=="X",23,chrom)) %>%
    mutate(chrom=as.integer(chrom))
  wbexp <- t(data.matrix(select(wbdf,-one_of(annocols))))
  write_df_h5(wbanno,group = "EXPinfo",outfile = h5filename,deflate_level = 4)

  write_2dmat_h5(h5f = h5filename,groupn = "EXPdata",datan = "expression",data = wbexp,deflate_level = 4,append=T)

  tf <- h5file(h5filename,'r')

  sexp <- tf['/EXPdata/expression'][,1:10]
  stopifnot(sexp==wbexp[,1:10])
  h5close(tf)
  return(TRUE)
}



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
  annodf <- select(all_genodata,one_of(annocols),doFlip) %>%
    select(-ref,-alt,-b37) %>%
    mutate(chrom=ifelse(chrom=="X",23,chrom)) %>%
    mutate(chrom=as.integer(chrom),pos=as.integer(pos))

  write_df_h5(df = annodf,group = "SNPinfo",outfile = h5file)
  rm(annodf)
  gc()
  all_genodata <-data.matrix(select(all_genodata,-one_of(annocols),-doFlip))
  gc()
  write_dmatrix_h5(h5file = h5file,groupname="SNPdata",dataname = "genotype",Nsnps = nrow(all_genodata),Nind = ncol(all_genodata),data = all_genodata,deflate_level = 4)
  rm(all_genodata)
  return(T)
}


