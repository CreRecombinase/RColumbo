
read_1kg_gz <- function(data_gzfile,leg_gzfile,map_gzfile,h5filename,subset_h5=NULL){
  require(dplyr)
  require(readr)
  require(h5)
  require(iotools)
  # data_gzfile <- "/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr5_1kg_geno.hap.gz"
  # leg_gzfile <- "/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr5_1kg_geno.legend.gz"
  # map_gzfile <- "/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/interpolated_from_hapmap/chr5.interpolated_genetic_map.gz"
  # h5filename <- "/media/nwknoblauch/Data/GTEx/1kg_SNP_H5/EUR.chr5_1kg.h5"
  # subset_h5 <- "/media/nwknoblauch/Data/GTEx/GTEx_SNP_h5/Adipose_Subcutaneous.h5"
  chr <- as.integer(gsub(".+EUR.chr([0-9]+)_1kg.+","\\1",data_gzfile))
  leg_data <- readr::read_delim(leg_gzfile,delim=" ",col_names = c("rsid","pos","ref","alt"),skip = 1) %>%
    mutate(snp_id=1:n(),chrom=chr)
  map_data <- readr::read_delim(map_gzfile,delim=" ",col_names = c("rsid","pos","map")) %>% mutate(chrom=chr,map_id=1:n())
  leg_map <- inner_join(leg_data,map_data) %>% distinct(snp_id,.keep_all=T) %>% distinct(map_id,.keep_all=T) %>% distinct(chrom,pos,.keep_all=T)
  w_leg_map <- select(leg_map,rsid,chrom,pos,ref,alt,map,snp_id)
  if(!is.null(subset_h5)){
    subset_leg <-read_df_h5(subset_h5,"SNPinfo") %>% filter(chrom==chr)
    w_leg_map <-  inner_join(w_leg_map,subset_leg,by=c("pos","chrom","ref","alt"))
  }
  w_leg_map <- distinct(w_leg_map,snp_id,.keep_all = T) %>% arrange(snp_id)
  write_df_h5(df = w_leg_map,outfile = h5filename,groupname = "SNPinfo",deflate_level = 4)
  write_rows_gzfile_h5(data_gzfile,w_leg_map$snp_id,h5filename,"SNPdata","genotype",max.size=335544320L)
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

  all_genodata <- data.table::fread(input = paste0("zcat ",ngzfile),sep = "\t",header = T,data.table = F)
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


