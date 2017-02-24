library(RColumbo)


args  <- commandArgs(trailingOnly = T)
gzfile <- args[1]
GDSfile <- args[2]
tissue <- gsub(".+/([^/]+)_Analysis.snps.txt.gz","\\1",gzfile)
gtissue <- gsub(".+/([^/]+).GDS","\\1",GDSfile)
stopifnot(tissue==gtissue)

write_gtex_gdsn(geno_txt = gzfile,geno_gdsn = GDSfile,chunksize = 650000)
