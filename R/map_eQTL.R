#Code for mapping eQTL

map_eqtl_lm <- function(Genotype,Expression,Covariates=NULL,scale_ortho_exp=F){
  require(dplyr)
  if(!is.null(Covariates)){
    covariates <- orthogonalize_covar(cbind(1,Covariates))
    ortho_Genotype <-orthogonalize_data(Genotype,covariates)
    if(is.null(dim(Expression))){
      ortho_Expression <- scale(orthogonalize_data(t(t(Expression)),covariates),scale = scale_ortho_exp,center = scale_ortho_exp)
    }else{
      ortho_Expression <- scale(orthogonalize_data(Expression,covariates),scale = scale_ortho_exp,center = scale_ortho_exp)
    }

  }else{
    ortho_Genotype=Genotype
    ortho_Expression=Expression
  }
#  summary(lm(ortho_Expression~ortho_Genotype[,1]+0))
  mmat <- eqtl_lm(ortho_Genotype,ortho_Expression)

  return(data_frame(betahat=mmat[,1],se=mmat[,2]))
}
