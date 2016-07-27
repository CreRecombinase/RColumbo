library(RColumbo)
context("LD Matrix Generation")

test_that("Reference haplotype data is correctly subsetted",{
rslistf <- "../rslist_test.txt"
legendfile <- "../test_1kg_geno_chr19.legend"
hapfile <- "../test_1kg_geno_chr19.hap"
mapfile <- "../test_OMNI_chr19_genetic_map.gz"
positionlistf= "../poslist_test.txt"
resfile_rs <- "../subset_test_rs.RDS"
resfile_pos <- "../subset_test_pos.RDS"
expect_equal(subset_ref_panel(rslistf = rslistf,
                              legendfile = legendfile,
                              mapfile = mapfile,
                              hapfile = hapfile),readRDS(resfile_rs))
expect_equal(subset_ref_panel(positionlistf = positionlistf,
                                          legendfile = legendfile,
                                          mapfile = mapfile,
                                          hapfile = hapfile),readRDS(resfile_pos))
})

test_that("subcols and subrows create subsets the way we want them to",{
  expect_equal(subcols(j = 1,chunksize = 4,ncols = 10),1:4)
  expect_equal(subrows(i=1,chunksize=4,ncols=10),1:4)
  expect_equal(subcols(j=1,chunksize=11,ncols=10),1:10)
  expect_equal(subcols(j=4,chunksize=4,ncols=15),13:15)
})
#
# test_that("LD matrices are equal to MATLAB implementation",{
#   reslist <-subset_ref_panel(rslistf = rslistf,
#                              legendfile = legendfile,
#                              mapfile = mapfile,
#                              hapfile = hapfile)
#   H <- t(reslist[["H"]])
#   cummap <- reslist[["cummap"]]
#   m <- 85
#   Ne <- 11490.672741
#   cutoff <- 1e-3
#   matresf <- "../true_R.csv"
#
#   true.R <- data.matrix(read.csv(matresf,header=F))
#
#
# })
# # test_that("LD generation works in parallel as well as in serial",{
  frsidf <- "~/Desktop/LDmapgen/fsnplist.txt"
  legendfile <- "~/Desktop/LDmapgen/tempdata/bmi2015_chr19_rss_1kg_geno_impute.impute.legend"
  hapfile <- "~/Desktop/LDmapgen/tempdata/bmi2015_chr19_rss_1kg_geno_impute.impute.hap"
  mapfile <- "~/Desktop/LDmapgen/1000-genomes-genetic-maps/interpolated_OMNI/chr19.OMNI.interpolated_genetic_map.gz"
  haph5 <- "~/Desktop/LDmapgen/1kgenotypes/IMPUTE/EUR.chr19_1kg_geno_hap.h5"
  frsid <- scan(frsidf,what=character())
   chr19 <- subset_ref_panel(rsids = frsid,legendfile = legendfile,
                             hapfile=hapfile,
                             mapfile = mapfile,outhapfile = haph5)
     tH <- chr19[["H"]]
     tmap <- chr19[["cummap"]]
    m <- 85
#
    cutoff <- 1e-3
    Ne <- 11490.672741
#   #
#   # # theta <
#
#   #
system.time(vecLD <- gen_LD(chr19,m = m,Ne=Ne,cutoff=cutoff))
ncor <- cov2cor(vecLD)
cLD <- as.matrix(nLD)
cLD <- cLD+t(cLD)
diag(cLD) <- 1
summary(c(cLD))
system.time(nLD <- p_sparse_LD(tmap, tH,  Ne,  m, cutoff,100))
system.time(nLD <- arm_gen_LD(tH,tmap,m,Ne,cutoff,5000))

  #
  #
  # sum(abs(nLD[upper.tri(nLD,T)]-corLD[upper.tri(corLD,diag = T)]))





  # system.time(sparse_LD
#   system.time(fLD <- fast_LD(tH,tmap,m,Ne,cutoff,4))
#
# })
