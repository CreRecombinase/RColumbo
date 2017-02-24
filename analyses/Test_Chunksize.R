
library(dplyr)
testfunction <- function(filechunksize,matrix,foth5file,nsnps){
  if(file.exists(foth5file)){
    file.remove(foth5file)
  }
  timing <- system.time(Rwrite_cov_LD(th5file = foth5file,tdimension = nsnps,Offset = 0,data = matrix,chunkvec = c(filechunksize,filechunksize)))[3]
  return(c(timing,file.info(foth5file)$size))
}

testreadfunction <- function(filechunksize,matrix,foth5file,nsnps){
  if(file.exists(foth5file)){
    file.remove(foth5file)
  }
  writetiming <- system.time(Rwrite_cov_LD(th5file = foth5file,tdimension = nsnps,Offset = 0,data = matrix,chunkvec = c(filechunksize,filechunksize)))[3]
  filechunksize <- min(2*filechunksize,ncol(matrix))
  readchunks = ceiling(ncol(matrix)/filechunksize)
  readtiming <-0
  for(i in 1:readchunks){
    for (j  in i:readchunks) {
      cat("i:",i,"\n")
      cat("j:",j,"\n")
      rowoffset <-(i-1)*filechunksize
      coloffset <- (j-1)*filechunksize
      treadtiming <- system.time(read_2dfmat_h5(h5file = foth5file,
                                                groupname = "LD_mat",
                                                dataname = "LD",
                                                row_offset = rowoffset,
                                                col_offset = coloffset,row_chunksize = filechunksize,col_chunksize = filechunksize))[3]
      readtiming=readtiming+treadtiming
    }
  }
  return(data_frame(readtiming=readtiming,writetiming=writetiming,filesize=file.info(foth5file)$size,chunksize=filechunksize))
}

chunksizes <- as.integer(seq(1001,12500,length.out = 10))
ntimings <- lapply(chunksizes,testreadfunction,matrix=fLDmat,foth5file=foth5file,nsnps=nsnps)
timedf <- bind_rows(timings)
plot(x=timedf$chunksize,y=timedf$readtiming)
plot(x=chunksizes,y=timings[1,])
plot(x=chunksizes,y=timings[2,])
