library(RColumbo)
library(Matrix)
library(dplyr)
library(readr)
library(rhdf5)
# library(future.BatchJobs)
# rsdir <- "~/Desktop/LDmapgen/rslists/"
# rsfiles <- dir(rsdir,full.names = T)
# rslen <- lapply(rsfiles,scan,what=character())
# rownum <- lengths(rslen)
chunksize <- 300000
#chunknums <- ceiling(rownum/chunksize)
#chunktot <- chunknums+((chunknums*chunknums)-chunknums)/2
mapfiles <- paste0("/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/interpolated_OMNI/chr",1:22,".OMNI.interpolated_genetic_map.gz")
hapfiles <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno.hap.gz")
legfiles <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno.legend.gz")
h5files <- paste0("/media/nwknoblauch/Data/1kg/haplotypes/EUR.chr",1:22,"_1kg_geno_hap.h5")
eqtldatfile <- paste0("~/Desktop/eQTL/Snake/Whole_Blood_Analysis.snpdata.txt.gz")
eqtllegfile <- paste0("~/Desktop/eQTL/Snake/WB_SNPS.txt")

for(i in 1:22){
  cat(paste0(i," of ",22,"\n"))
  hapfile <- hapfiles[i]
  mapfile <- mapfiles[i]
  h5file <- h5files[i]
  legfile <- legfiles[i]
  legdf <-read.table(legfile,header=T,sep=" ",stringsAsFactors = F)
  nSNPs <- nrow(legdf)
  ncols <- length(scan(hapfile,what=integer(),sep=" ",nlines = 1))
  write_haplotype_h5(hapfile,h5file,nSNPs,ncols = ncols,chunksize = 250000,deflate_level = 3)
  write_leg_h5(legendfile = legfile,haph5 = h5file,chunksize = 125000)
  write_map_h5(mapfile=mapfile,haph5=h5file,chunksize=125000)
}

legendfile <- "~/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno.legend.gz"
hapfile <- "/home/nwknoblauch/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno.hap.gz"
mapfile <- "~/Desktop/LDmapgen/1000-genomes-genetic-maps/interpolated_OMNI/chr19.OMNI.interpolated_genetic_map.gz"
haph5 <- "/home/nwknoblauch/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno_hap_2.h5"
eqtlfile <- "~/Desktop/eQTL/Snake/IBD_WholeBlood_eQTL.h5"
result_dir <- "/home/nwknoblauch/Desktop/LDmapgen/"


mapdf <- read_h5_df(haph5,"Map")
rslistvec <- subset_h5_ref(haph5,mapfile = mapfile,eqtlfile)
nSNPs <- length(rslistvec)
chunksize <- 30000
chunknum <- ceiling(nSNPs/chunksize)
chrom <- 19
m=85
Ne=11490.672741
cutoff=1e-3
i <- 0
j <- 0
bsl <- list()
k <- 1
ct <- Sys.time()
nt <- Sys.time()
for(i in 0:(chunknum-1)){
  for(j in i:(chunknum-1)){
    bsn <-paste0(i,"_",j)
    if(is.null(bsl[[paste0(i,"_",j)]])){
      if(!file.exists(bsfile <- paste0("~/Desktop/LDmapgen/chr19/chr",chrom,"_",i,"_",j,"_",chunksize,".RDS"))){
        cat(paste0(i,"_",j,"_",nt-ct,"\n"))
        ct <- Sys.time()
        cat(paste0("\n",Sys.time()))
        bsl[[paste0(i,"_",j)]] <- torque_arm_gen_LD(rslistvec = rslistvec,eqtlfile = eqtlfile,
                                                    haph5 = haph5,
                                                    mapfile = mapfile,
                                                    m = m,Ne = Ne,cutoff = cutoff,
                                                    chunksize = chunksize,
                                                    i = i,j = j,
                                                    result_dir = result_dir,
                                                    chrom = chrom)
        nt <- Sys.time()
      }
    }
  }
}
mapdf <-read.table(mapfile,header=F,stringsAsFactors = F) %>% select(rsid=V1,pos=V2,cummap=V3)
mapdf <- filter(mapdf,rsid %in% rslistvec)
cummap <- mapdf$cummap
write.table(cummap,"~/Desktop/LDmapgen/bigtest_cummap.txt",col.names=F,row.names=F,sep=" ")
matlab_r <- readMATLAB("~/Desktop/LDmapgen/testRColumbo.h5")


bmat <- new("dgTMatrix",i=il,j=unlist(jl),x=unlist(xl),Dim=dim(mats[[1]]))






tdbs <- read_h5_df("~/Desktop/eQTL/Snake/dbsnp.h5","dbSNP")
head(tdbs)



 H5read_haplotype_h5(haph5,1:10)
 rsvec <- read_rsid_h5(haph5)
 rslist <-subset_h5_ref(haph5 = haph5,mapfile = mapfile,eqtlfile =eqtlfile)
 fmat <- p_arm_gen_LD(rsids,haph5 = haph5,mapfile = mapfile,m = m,Ne = Ne,cutoff = cutoff,chunksize = 10000)


 hl <-subset_ref_panel(rsids=rslist,legendfile=legendfile,
                       hapfile=hapfile,
                       mapfile=mapfile,outhapfile = haph5)

wbs <- read.table("/media/nwknoblauch/Data/GTEx/GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices/WB_SNPS.txt",header=F,skip = 1,sep="_",col.names=c("chr","pos","ref","alt","b"),stringsAsFactors = F)

