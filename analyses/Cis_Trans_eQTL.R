#Code to try to reconstruct eQTL
library(RColumbo)
library(rhdf5)
library(dplyr)
library(tidyr)
#eqtlh5 <- "~/Desktop/eQTL/Snake/IBD_WholeBlood_eQTL.h5"
rawh5 <- "/home/nwknoblauch/Desktop/eQTL/Snake/Whole_Blood_eQTL_raw_data.h5"
oh5s <- paste0("/media/nwknoblauch/Data/GTEx/chr",1:22,"_Whole_Blood_trans_eQTL.h5")
haph5s <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
gwash5 <- "~/Desktop/eQTL/Snake/IBD.h5"
vegasf <- "~/Downloads/VEGAS/EUR.IBD_filtered_Vegas_results.txt.out"
gmf <- "/media/nwknoblauch/Data/GTEx/gencode.v19.genes.patched_contigs.gtf.gz"
# snpleg <- read_h5_df(rawh5,"SNPinfo") %>% mutate(snp_ind=1:length(rsid))
expleg <- as_data_frame(t(h5read(rawh5,"EXPinfo/annomat"))) %>% rename(fgeneid=V1,sgeneid=V2,chrom=V3,TSStart=V4,TSStop=V5) %>% mutate(exp_id=1:(length(fgeneid)))

gc()
gmd <-vegas.merge(vegasf,gmf)
gmd <- filter(gmd,!duplicated(ensid))
gmd <- mutate(gmd,fgeneid=as.integer(gsub("ENSG","",ensid)))
expleg <- inner_join(expleg,gmd,by="fgeneid")
ilrleg <- filter(expleg,Gene=="IL23R")
subleg <- ilrleg



#sub_snpleg <- subset_all_hap(snpleg,haph5s,gwash5)
#saveRDS(sub_snpleg,"~/Desktop/eQTL/hap_eQTL.RDS")
sub_snpleg <- readRDS("~/Desktop/eQTL/hap_eQTL.RDS") %>% select(pos,rsid,chrom,hsnpid,cummap,snp_ind) %>% mutate(rsidi=as.integer(gsub("rs","",rsid)))

hone <- filter(sub_snpleg,chrom==1) %>% select(hsnpid)
hone <- hone[[1]]
dokeep <- as.integer(1:1891451 %in% hone)
mhmat <- read_fmat_gz("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr1_1kg_geno.hap.gz",chunksize = 750000,nrows = 1891451,keeprow = dokeep,ncols = 1006)
mhmat <- t(mhmat)
thmat <- read_haplotype_ind_h5(haph5s[1],hone)
oleg <- read_h5_df(haph5s[1],"Legend")
legs <- slice(oleg,hone)
sum(abs(mhmat-thmat))








eqtlh5 <- "~/Desktop/eQTL/Snake/IBD_WholeBlood_eQTL.h5"
h5ls(eqtlh5)
eqtl_df <- read_h5_df(eqtlh5,"eQTL")
eqtl_df <- select(eqtl_df,-sgeneid,-betahat,-pvalg,-seg)
snp_df <- read_h5_df(eqtlh5,"SNP")
nsnp_df <- select(sub_snpleg,-rsid) %>%rename(rsid=rsidi) %>% inner_join(snp_df)
eqtl_df <- inner_join(eqtl_df,nsnp_df)
head(eqtl_df)
ilr_eqtl <- inner_join(ilrleg,eqtl_df)
head(ilr_eqtl)
haph5 <- haph5s[1]
fmat <- comp_dense_LD(ilr_eqtl$hsnpid,ilr_eqtl$cummap,ilr_eqtl$rsid,haph5 = haph5,m=85,Ne=11490.672741,cutoff=1e-3)
system.time(tfmat <- gen_dense_LD(haph5,ilr_eqtl$hsnpid,ilr_eqtl$cummap,m=85,Ne=11490.672741,cutoff=1e-3,0,0,length(ilr_eqtl$hsnpid)))
system.time(ofmat <- pgen_dense_LD(haph5,ilr_eqtl$hsnpid,ilr_eqtl$cummap,m=85,Ne=11490.672741,cutoff=1e-3,0,0,length(ilr_eqtl$hsnpid)))
subgraphs <- find_subgraphs(fmat,0.5,ilr_eqtl$rsid)
tl <- list()
ti <- 1
# for(i in 1:length(subgraphs)){
#   for(j in 1:length(subgraphs)){
#     cat(paste0("i:",i,"j:",j,"\n"))
#     sfmat <- fmat[as.character(subgraphs[[i]]),as.character(subgraphs[[j]])]
#     tl[[ti]] <- data_frame(i=i,j=j,size=length(sfmat),nconnected=sum(abs(sfmat[sfmat!=1])>0.5),disconnected=all(abs(sfmat)<0.5))
#     ti <- ti+1
#   }
# }
tl <- list()
for(i in 1:length(subgraphs)){
tl[[i]] <- data_frame(subgraph=i,rsid=subgraphs[[i]])
}
adf <- bind_rows(tl)
nirl <- inner_join(ilr_eqtl,adf)
bestb <- group_by(nirl,subgraph) %>% filter(pvale==min(pvale))

lapply(subgraphs)
tldf <- bind_rows(tl)
filter(tldf,i!=j) %>% summarise(nd=sum(disconnected))



hist(c(abs(fmat)))





#bexpdat <- read_dmat_ind_h5(rawh5,"EXPdata","orthoexpression",c(1:23973))

my_db <- src_sqlite(path = "~/Desktop/eQTL/SNP_exp.db", create = F)
nsub_snpleg <- tbl(my_db,"SNP")
#nsub_snpleg <- copy_to(my_db,sub_snpleg,name="SNP",temporary=F,indexes=list(c("chrom","pos"),"pos"))
#nexp <- copy_to(my_db,expleg,"EXP",temporary=F,indexes = list(c("chrom","TSStart"),c("chrom","TSStop"),"TSStart","TSStop"))
nexp <- tbl(my_db,"EXP")
group_by(nexp,chrom) %>% summarise(ngenes=n())
# for(i in 1:22){
#   thaph5 <- haph5s[i]
#   oh5 <- oh5s[i]
#   tchsnp <- filter(nsub_snpleg,chrom==i)
#   tchexp <- filter(nexp,chrom==i)
#   ugenes <- collect(distinct(tchexp,fgeneid))[["fgeneid"]]
#   for(j in 2:length(ugenes)){
#     cat(paste0(j," of ",length(ugenes),"\n"))
#     ciseq <-inner_join(tchexp,tchsnp,by="chrom") %>% mutate(iscis=((pos<TSStop&pos>TSStart)|(abs(pos-TSStart)<1e6)|(abs(pos-TSStop)<1e6))+0)
#     ciseq <- arrange(ciseq,snp_ind,exp_id)
#     dfciseq <- collect(ciseq)
#
#     nsnps <- length(ciseq$fgeneid)
#     if(nsnps>0){
#     cis_stat_extract(query_df = ciseq,haph5 = thaph5,rawh5 = rawh5,outh5 = oh5,chromosome = i,pcutoff = 0.01,LDcutoff = 0.8,append = T,chunksize=10000 )
#     }
#   }
# }











outh5 <- "/media/nwknoblauch/Data/GTEx/chr19_Whole_Blood_trans_eQTL.h5"
tcutoff <- 3.5
LDcutoff <- 0.8
stat_extract(eqtlh5,rawh5,haph5,gwash5,outh5,10000,tcutoff,LDcutoff)




