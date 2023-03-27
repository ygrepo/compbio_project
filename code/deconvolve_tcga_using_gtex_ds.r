rm(list = ls())
# load
library(MuSiC)
library(SingleCellExperiment)
library(dplyr)
#library(nnls)

setwd("~/github/compbio_project/code")

# source("music/utils.R")
# source("music/construct.R")
# source("music/analysis.R")

tcga_file = paste("../data/brca/tcga/processed/", 
                  "converted_genes_primary_tumor_unstranded_exp_data.rds", 
                  sep="")
exprs <-readRDS(tcga_file)
#exprs <- exprs[,1:2]
print(dim(exprs))

sce_file = paste("../data/brca/gtex/processed", 
                 "subset_sce.rds", 
                 sep="")
sce <- readRDS(sce_file)
print(dim(sce))

ct <- unique(sce$CellType)
print(ct)

sce_genes <- toupper(rowData(sce)$value)
tcga_genes <- toupper(rownames(exprs))
rownames(exprs) <- tcga_genes

markers <- intersect(sce_genes, tcga_genes)
msg <- paste("Markers", length(markers), sep=":")
print(msg)
msg <- paste("#sce genes", length(sce_genes), sep=":")
print(msg)
msg <- paste("#tcga genes", length(tcga_genes), sep=":")
print(msg)


# debug(music_prop)
# Estimate cell type proportions
prop = music_prop(bulk.mtx = exprs, 
                                markers = markers,
                                sc.sce = sce, 
                                clusters = 'CellType', 
                                samples = 'sampleID',
                                select.ct = ct,
                                verbose = TRUE)
prop_file = paste("../data/brca/tcga/processed/gtex/", 
                  "prop_primary_tumor_unstranded.rds", 
                  sep="")
msg = paste("Saving props to ", prop_file, sep="")
print(msg)
saveRDS(prop, file=prop_file)

#undebug(music_prop)