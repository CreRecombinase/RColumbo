                                        #Generate LD matrices
library(methods)
library(RColumbo)
library(rhdf5)
library(dplyr)
                                        # library(ggplot2)
chrom <- 22
chrom <- as.integer(commandArgs(trailingOnly=T))
h5files <- paste0("/group/xhe-lab/RSS/EUR.chr",1:22,"_1kg_Whole_Blood_HEIGHT.h5")
#th5file <- "/group/xhe-lab/RSS//media/nwknoblauch/Data/GTEx/chr22_1kg_Whole_Blood_HEIGHT.h5"
oh5files <- paste0("/group/xhe-lab/RSS/LD_matrices/chr",1:22,"_1kg_WB_Height_LD.h5")
foh5files <-oh5files

# get_rownum_h5(h5files[2],"Haplotype","genotype")
# chunksize <- get_rownum_h5(hap_h5file = th5file,"Haplotype","genotype")
## chrom <- 22
 chunksize <- 1250
## cat(chrom,"\n")
## th5file <- h5files[chrom]
## foth5file <- foh5files[chrom]
## # oth5file <- oh5files[chrom]
## nsnps <-get_rownum_h5(hap_h5file = th5file,"Haplotype","genotype")
## chunknum <- ceiling(nsnps/chunksize)
## offset <-0
##  cLD <- read_2dfmat_h5(foth5file,"LD_mat","LD",row_offset = 0,col_offset = 0,row_chunksize = chunksize*2,col_chunksize = chunksize*2)
## cLD <- read_2dfmat_h5(foth5file,"LD_mat","LD",row_offset = 0,col_offset = 25000-1,row_chunksize = chunksize,col_chunksize = chunksize)

#for(chrom in 1:21){
cat(chrom,"\n")
th5file <- h5files[chrom]
foth5file <- foh5files[chrom]
stopifnot(file.exists(th5file))
                                        # oth5file <- oh5files[chrom]
nsnps <-get_rownum_h5(hap_h5file = th5file,"Haplotype","genotype")
chunknum <- ceiling(nsnps/chunksize)
offset <-0

#thapdat <- read_fmat_h5(th5file,"Haplotype","genotype",0,500)

if(file.exists(foth5file)){
    file.remove(foth5file)
}


testfunction <- function(chunksize,matrix){
    

for(j in 1:chunknum){
    cat(j," of ",chunknum,"\n")
    fLDmat <- fslide_LD(th5file,chunksize,offset,1e-3)
Rwrite_blosc_cov_LD(th5file = foth5file,tdimension = nsnps,Offset = offset,data = fLDmat)
    offset <- offset+chunksize
}
                                        #}

## cLD <- read_2dfmat_h5(foth5file,"LD_mat","LD",row_offset = 0,col_offset = 0,row_chunksize = chunksize*2,col_chunksize = chunksize*2)

## # system.time(ccLD <- read_2dfmat_h5(foth5file,"LD_mat","LD",0,0,row_chunksize = nsnps,col_chunksize = nsnps))
## # ccLD[1:5,1:5]
## # cLD[1:5,1:5]
## # head(which(fLDmat[1:nrow(cLD),1:ncol(cLD)]!=cLD,arr.ind=dim(cLD)))
## # head(which(ccLDmat[1:nrow(fLDmat),1:ncol(fLDmat)]!=fLDmat,arr.ind = dim(fLDmat)))
## # sum(fLDmat[(nrow(ofLDmat)+1):nrow(fLDmat),(nrow(ofLDmat)+1):nrow(fLDmat)]!=tfLDmat)/length(tfLDmat)

## # fLDmat <- fslide_LD(th5file,chunksize,0,1e-3)
## tadjmat <- find_adjmat_chunk(foth5file,0.25,12500)
## chrom <- 22
## mapdf <- read_h5_df(h5files[chrom],"Legend")
## ccLDmat <- h5read(oth5file,name = "/LD_mat/LD")
## ccLDmat <- ccLDmat+t(ccLDmat)
## indLDmat <- find_adjmat(ccLDmat,0.25,0,0)
## tindLDmat <- find_adjmat(fLDmat,0.25,0,0)
## colormat <- sequential_coloring(tindLDmat)
## colordf <- data_frame(ind=colormat[,1]+1,color=colormat[,2])
## subind <- filter(colordf)
## ggplot(colordf)+geom_histogram(aes(x=color))

## plot(ccLDmat[3,])
## head(tfLDmat[29060-25000+1,])
## head(ccLDmat[29059,])-head(fLDmat[29059,])
## subind <- find_adjmat(subLD,0.25,0,10000)



## # cLDmat[1:5,1:5]
## # fLDmat[1:5,1:5]
## # sum(abs(fLDmat-cLDmat))



