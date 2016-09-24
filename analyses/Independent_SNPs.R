library(RColumbo)
library(rhdf5)
library(dplyr)
library(tidyr)
oh5s <- paste0("/media/nwknoblauch/Data/GTEx/chr",1:22,"_1kg_Whole_Blood_IBD.h5")
for(chrom in 3:22){
  cat(paste0("chr",chrom,"\n"))
  h5file <- oh5s[chrom]
  map <- c(h5read(h5file,"/Legend/cummap"))
  m=85
  Ne=11490.672741
  cutoff=1e-3
  chunksize <- 10000

# system.time(mtfmat <- chunk_LD(h5file = h5file,tmap = map[1:30000],offset = 0,chunksize = 30000,m = m,Ne = Ne,cutoff = cutoff))
# sptfmat <- slide_LD(h5file = h5file,tmap = map,offset = 0,chunksize = chunksize,m = m,Ne = Ne,cutoff = cutoff)
# nsptfmat <-slide_LD(h5file = h5file,tmap = map,offset = chunksize,chunksize = chunksize,m = m,Ne = Ne,cutoff = cutoff)
# sum(abs(mtfmat[1:10000,1:20000]-sptfmat))
# sum(abs(mtfmat[10001:20000,10001:30000]-nsptfmat))
# ifmata <- cutoff_LD(sptfmat,offset = 0,LD_cutoff = 0.25)
# ifmatb <- cutoff_LD(nsptfmat,offset = chunksize,LD_cutoff = 0.25)

  offset <- 0
  tfml <- list()
  j <- 1
  cat(paste0(j))
  tfml[[j]]<- slide_cutoff_LD(h5file = h5file,tmap = map,offset = offset,chunksize = chunksize,m = m,Ne = Ne,cutoff = cutoff,LD_cutoff = 0.25)
  #compl[[j]]<-find_subgraphs(tedges = tfml[[j]],n_vertices =chunksize*2,offset = offset)
  offset <- offset+chunksize
  j <- j+1


  while(max(tfml[[j-1]][,1])+1<length(map)){
    cat(paste0(j,"\n"))
    tfml[[j]]<- slide_cutoff_LD(h5file = h5file,tmap = map,offset = offset,chunksize = chunksize,m = m,Ne = Ne,cutoff = cutoff,LD_cutoff = 0.25)
    #  compl[[j]]<-find_subgraphs(tedges = tfml[[j]],length(unique(c(tfml[[j]]))))
    offset <- offset+chunksize
    j <- j+1
  }
  apairs <- do.call("rbind",tfml)
  n_elem <- length(unique(apairs[,1]))
  saveRDS(apairs,paste0("/media/nwknoblauch/Data/subsets_chr",chrom,".RDS"))
  comp <- find_subgraphs(apairs,n_vertices = n_elem,offset = 0)
  compdf <- data_frame(g_ind=comp[,1],group=comp[,2])
  gwas_leg <- read_h5_df(h5file,"Legend") %>% mutate(g_ind=1:n())
  ngwas_leg <- inner_join(gwas_leg,compdf)
  saveRDS(ngwas_leg,paste0("/media/nwknoblauch/Data/IBD_GWAS_chr",chrom,".RDS"))
}

