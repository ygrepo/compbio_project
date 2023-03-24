rm(list = ls())
# load
# library(AnnotationDbi)
# library(EnsDb.Hsapiens.v79)
#library(MuSiC)
library(SingleCellExperiment)
library(dplyr)
#library(nnls)

setwd("~/github/compbio_project/code")

# source("music/utils.R")
# source("music/construct.R")
# source("music/analysis.R")

tcga_file = paste("../data/brca/tcga/processed/", 
                  "converted_genes_normal_tissue_unstranded_exp_data.rds", 
                  sep="")
exprs <-readRDS(tcga_file)


sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "sce_major_cell_type.rds", 
                 sep="")
sce <- readRDS(sce_file)
ct <- unique(sce$cellType)



#debug(music_prop)
# Estimate cell type proportions
prop = music_prop(bulk.mtx = exprs, 
                                sc.sce = sce, 
                                clusters = 'cellType', 
                                samples = 'Barcode',
                                select.ct = ct,
                                verbose = TRUE)
prop_file = paste("../data/brca/tcga/processed/Wu/", 
                        "prop_normal_tissue_unstranded_major_ct.rds", sep="")
saveRDS(prop, file=prop_file)

#undebug(music_prop)