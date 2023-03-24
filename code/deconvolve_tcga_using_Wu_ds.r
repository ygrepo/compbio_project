rm(list = ls())
# load
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
#library(MuSiC)
library(SingleCellExperiment)
library(dplyr)
library(nnls)

setwd("~/github/compbio_project/code")

source("music/utils.R")
source("music/construct.R")
source("music/analysis.R")

tcga_file = paste("../data/brca/tcga/processed/", 
                  "convert_genes_unstranded.rds", sep="")
tcga <-readRDS(tcga_file)


exprs <- as.matrix(tcga[, -1])
# take only two first samples
exprs <- exprs[,1:2]
rownames(exprs) <- rownames(tcga)


sc_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "small_wu_etal_sce.rds", sep="")
sc_data <- readRDS(file=sc_file)
ct <- unique(sc_data$cellType)



#debug(music_prop)
# Estimate cell type proportions
# debug(music_prop)
prop = music_prop(bulk.mtx = exprs, sc.sce = sc_data, 
                                clusters = 'cellType', 
                                samples = 'Barcode',
                                select.ct = ct,
                                verbose = TRUE)
prop_file = paste("../data/brca/tcga/processed/Wu/", 
                        "03_24_small_tcga_small_sc_prop.csv", sep="")
saveRDS(prop_cells, file=Est.prop.TCGA_prop)

#undebug(music_prop)