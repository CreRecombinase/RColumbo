library(sparklyr)
library(dplyr)
library(RColumbo)

sc <- spark_connect(master = "local")

gtex_h5 <- "~/Desktop/eQTL/Snake/Whole_Blood_eQTL_v6p_raw_data.h5"
idgtexf <- "/media/nwknoblauch/Data/GTEx/ind_gtex_snps.RDS"
ind_gtex_leg <- readRDS(idgtexf) %>% rename(pos=position)
gwas_tbl<- copy_to(sc,gwas_dat)

gwasRDS <- "/media/nwknoblauch/Data/gwas_data/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc.RDS"
gwas_dat <- readRDS(gwasRDS)
gwas_dat <- mutate(gwas_dat,doFlip=as.integer(A1>A2),
                   nOR=ifelse(doFlip==1,1/OR,OR),
                   posl=strsplit(Direction,split=""))%>%
  mutate(npos=map_int(posl,~sum(.=="+")),
         nneg=map_int(posl,~sum(.=="-")),
         namb=map_int(posl,~sum(.=="?"))) %>% dplyr::select(-posl)
