library(RColumbo)
library(rhdf5)
gwdh <- strsplit(getwd(),split = "/",fixed=T)[[1]][3]
if(gwdh=="nwknoblauch"){
   oh5files <- paste0("/media/nwknoblauch/Data/DGN/LDmats/chr",1:22,"_1kg_DGN_Height_LD.h5")
  # oh5files <- paste0("/media/nwknoblauch/Data/GTEx/LDmats/chr",1:22,"_1kg_WB_Height_LD.h5")
}else{
  h5files <- paste0("/group/xhe-lab/RSS/EUR.chr",1:22,"_1kg_Whole_Blood_HEIGHT.h5")
  oh5files <- paste0("/group/xhe-lab/RSS/LD_matrices/chr",1:22,"_1kg_WB_Height_LD.h5")
}

chrom <- 1
for(chrom in 2:22){
  # chrom <- as.integer(commandArgs(trailingOnly = T))
  cat(chrom,"\n")
  moh5file <- oh5files[chrom]
  adjdat <- rbind(read_Rint_h5(moh5file,"LD_adj","rowid"),
                  read_Rint_h5(moh5file,"LD_adj","colid"))
  sum(adjdat[1,]!=adjdat[2,])
  color <- sequential_coloring(adjdat)
  write_Rint_h5(moh5file,"Color","index",color[,1],4)
  write_Rint_h5(moh5file,"Color","color",color[,2],4)
}

