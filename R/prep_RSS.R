
# param_snpfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Whole_Blood_1KG/Adipose_Subcutaneous/SNP.chunks.h5"
# snp_chunk <- 1
# exp_chunk <-1
# param_expfile <- "/media/nwknoblauch/Data/GTEx/GTEx_Whole_Blood_1KG/Adipose_Subcutaneous/EXP_Adipose_Subcutaneous.h5"
# tissue <- "Adipose_Subcutaneous"
# scale_ortho_exp=T
# m=85
# Ne=11490.672741
# cutoff=1e-3
# em_logodds=F
run_RSS <- function(param_snpfile,param_expfile,snp_chunk,exp_chunk,tissue,
                    scale_ortho_exp=F,m=85,Ne=11490.672741,cutoff=1e-3,
                    logodds=NULL,sigb=NULL,exp_chunksize=1){
  require(Matrix)
  require(dplyr)
  require(rssr)
  library(h5)
  library(BBmisc)
  stopifnot(file.exists(param_snpfile),
            file.exists(param_expfile))


  eqtl_snpfile <-read_attr(param_snpfile,as.character(snp_chunk),paste0(tissue,"_filepath"))
  eqtl_snp_chromosome <-read_attr(param_snpfile,as.character(snp_chunk),"chromosome")
  stopifnot(file.exists(eqtl_snpfile))

  eqtl_snplist <-read_vec(param_snpfile,paste0(snp_chunk,"/",tissue))
  eqtl_snpdat <- read_ind_h5(eqtl_snpfile,"SNPdata","genotype",eqtl_snplist)



  if(exp_chunksize==1){
    eqtl_explist <-read_vec(param_expfile,paste0("all/",tissue))[exp_chunk]
    exp_fgeneid <- read_vec(param_expfile,"all/fgeneid")[exp_chunk]
    eqtl_expfile <- read_attr(param_expfile,"all",paste0(tissue,"_filepath"))
  }else{
    if(!exists_group(param_expfile,as.character(exp_chunk))){
      eqtl_explist <-chunk(read_vec(param_expfile,paste0("all/",tissue)),chunk.size = exp_chunksize)[[exp_chunk]]
      eqtl_fgeneid <-chunk(read_vec(param_expfile,paste0("all/fgeneid")),chunk.size = exp_chunksize)[[exp_chunk]]
      eqtl_expfile <- read_attr(param_expfile,"all",paste0(tissue,"_filepath"))
    }else{
    eqtl_explist <-read_vec(param_expfile,paste0(exp_chunk,"/",tissue))
    exp_fgeneid <- read_vec(param_expfile,paste0(exp_chunk,"/fgeneid"))
    eqtl_expfile <- read_attr(param_expfile,as.character(exp_chunk),paste0(tissue,"_filepath"))
    }
  }
  stopifnot(file.exists(eqtl_expfile))
  eqtl_expdat <- read_ind_h5(eqtl_expfile,"EXPdata","expression",eqtl_explist)


  eqtl_covdat <- read_2d_mat_h5(eqtl_expfile,"Covardat","covariates")

  teqtl_df <- list()
  for(i in 1:ncol(eqtl_expdat)){
    teqtl_df[[i]] <-map_eqtl_lm(eqtl_snpdat,eqtl_expdat[,i],eqtl_covdat,scale_ortho_exp = scale_ortho_exp) %>% mutate(fgeneid=exp_fgeneid[i])
  }
  summary_stats <- bind_rows(teqtl_df)

  ld_ldf<-read_attr(param_snpfile,as.character(snp_chunk),"1kg_filepath")
  stopifnot(file.exists(ld_ldf))
  ld_snplist <- read_vec(param_snpfile,paste0("/",snp_chunk,"/1kg"))

  snpA <-read_ind_h5(ld_ldf,"SNPdata","genotype",ld_snplist)
  mapA <-read_vec(ld_ldf,"/SNPinfo/map")[ld_snplist]

  stopifnot(ncol(snpA)==length(ld_snplist),length(mapA)==length(ld_snplist))


  R <- gen_sparsemat(calcLD(hmata=snpA,hmatb = snpA,mapa = mapA,mapb = mapA,m = m,Ne=Ne,cutoff = cutoff,isDiag = T),istart=1,jstart=1,nSNPs = ncol(snpA),makeSymmetric = T)
  rm(snpA,mapA,eqtl_snpdat,eqtl_expdat,eqtl_covdat)

  fit_dfl <- list()
  summ_statl <- split(summary_stats,summary_stats$fgeneid)
  for(i in 1:length(summ_statl)){
    cat(i,"of ",length(summ_statl),"\n")
    se <- summ_statl[[i]][["se"]]
    betahat <- summ_statl[[i]][["betahat"]]
    SiRiS <- SiRSi(R,Si = 1/se)
    p <- length(betahat)
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- SiRiS%*%(alpha0*mu0)
    fit_df <- grid_search_rss_varbvsr(SiRiS = SiRiS,sigma_beta = sigb,logodds=logodds,betahat = betahat,
                                                                      se = se,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = SiRiSr@x,tolerance = 1e-3,itermax = 100,verbose = T,lnz_tol = T)

    fit_dfl[[i]] <- fit_df %>% mutate(fgeneid=exp_fgeneid[i])
    gc()
  }
  fit_df <- bind_rows(fit_dfl) %>% mutate(snp_chunk=snp_chunk,exp_chunk=exp_chunk,snp_chromosome=eqtl_snp_chromosome)
  return(fit_df)
}
