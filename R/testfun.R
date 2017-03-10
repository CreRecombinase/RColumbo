#basic implementations of RColumbo functions for testing purposes


test_eqtl <- function(X,Y,scale_x=F,scale_y=T){
  n <- nrow(X)
  X <- scale(X,center = scale_x,scale = scale_x)
  Y <- scale(Y,center=scale_y,scale=scale_y)
  res_array <-array(data=NA,dim=c(ncol(X),ncol(Y),2))

  for(i in 1:ncol(X)){
    for(j in 1:ncol(Y)){
      tlm <- coef(summary(lm(Y[,j]~X[,i]+0)))
      res_array[i,j,1] <- tlm[1]
      res_array[i,j,2] <- tlm[2]
    }
  }
  return(res_array)
}

test_betahat <- function(X,Y){
  res <- test_eqtl(X,Y)[,,1]
  return(res)
}

test_betahat <- function(x,y){
  xtx <- x%*%x
  n <- length(x)
  k <- 1
  b <-(x%*%y)/xtx
  r <- y-x*b
  s2 <- (r%*%r)/(n-k)
  serr <- sqrt(s2/xtx)

}

# sxyf <- function(X,Y){
#
# }

# test_eqtl_p <- function(X,Y,scale_x=F,scale_y=T){
#
#   n <- nrow(X)
#   X <- scale(X,center = scale_x,scale = scale_x)
#   Y <- scale(Y,center=scale_y,scale=scale_y)
#
#   sxx <-apply(X,2,function(x)sum((x-mean(x))^2))
#   syy <- apply(Y,2,function(x)sum((x-mean(x))^2))
#   sxy <- crossprod(scale(X),Y)
#   b <- sxy/sxx
#   se <- outer(syy,sxy,function(y,x){
#     sqrt(((y-x)/(n-2))/x)
#   })
#   se <- sqrt(((syy-sxy)/(n-2))/sxx)
#   res_array <-array(data=NA,dim=c(ncol(X),ncol(Y),2))
#   res_array[,,0] <- b
#   res_array[,,1] <- se
#   return(res_array)
# }

test_cov <- function(A,B,isDiag=T){
  if(isDiag){
    retmat <- cov(A,B)
    retmat[lower.tri(retmat)] <- 0
    return(retmat)
  }else{
    return(cov(A,B))
  }
}

test_dist <- function(mapa,mapb, isDiag){

  distmat <- matrix(0,nrow = length(mapa),ncol = length(mapb))
  for(i in 1:nrow(distmat)){
    for(j in min(c(ncol(distmat),(i+1))):ncol(distmat)){
      distmat[i,j] <- mapb[j]-mapa[i]
    }
  }
  return(distmat)
}

# reqtlf <-'/media/nwknoblauch/Data/GTEx/GTEx_cis_eQTL_h5/Whole_Blood.h5'
