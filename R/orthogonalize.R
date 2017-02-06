
orthogonalize_dataset <- function(h5filename,newh5filename,covar_h5file,datagroup,datasetname,newdatasetname,chunksize){

  require(h5)
  require(BBmisc)
  cvdat <- read_2d_mat_h5(covar_h5file,groupn = "Covardat","covariates")
  cvdat <- cbind(1,cvdat)
  covariates <- orthogonalize_covar(cvdat)

  dataf <- h5file(h5filename,mode='r')
  datag <- dataf[datagroup]
  datad <- datag[datasetname]
  nitems <- ncol(datad)
  chunkseq <- chunk(1:nitems,chunk.size = chunksize)
  nchunks <- length(chunkseq)

  odataf <- h5file(newh5filename,mode='a')
  odatag <- odataf[datagroup]
  odatad <- createDataSet(odatag,datasetname = newdatasetname,type = "double",
                          dimensions = as.integer(dim(datad)),
                          chunksize = as.integer(c(nrow(datad),5000L)),maxdimensions = as.integer(dim(datad)),compression = 4L)
  for(i in 1:nchunks){
    cat("Chunk ",i," of ",nchunks,"\n")
    odatad[,chunkseq[[i]]] <- orthogonalize_data(datad[,chunkseq[[i]]],covariates)
  }
  h5close(dataf)
  h5close(odataf)
}

