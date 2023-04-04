rm(list = ls())
# load
# library(AnnotationDbi)
# library(EnsDb.Hsapiens.v79)
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
msg <- paste("Tcga exp. dimensions", dim(exprs), sep=":")
print(msg)

sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/ct_minor/", 
                 "subset_CID3586.rds", 
                 sep="")
msg <- paste("Loading",
             sce_file,
             sep=" ")
print(msg)
sce <- readRDS(sce_file)
msg <- paste("Sce dimensions", dim(sce), sep=":")
print(msg)

#sce <- sce[,1:1000]
ct <- unique(sce$cellType)
msg <- paste("Cell type", ct, sep=":")
print(msg)

tcga_genes <- toupper(rownames(exprs))
rownames(exprs) <- tcga_genes

sce_genes <- toupper(rownames(sce))
rownames(sce) <- sce_genes

markers <- intersect(sce_genes, tcga_genes)
msg <- paste("Markers", length(markers), sep=":")
print(msg)
msg <- paste("#Sce genes", length(sce_genes), sep=":")
print(msg)
msg <- paste("#Tcga genes", length(tcga_genes), sep=":")
print(msg)


#debug(music_prop)
# Estimate cell type proportions
prop = music_prop(bulk.mtx = exprs,
                                markers = markers,
                                sc.sce = sce, 
                                clusters = 'cellType', 
                                samples = 'Barcode',
                                select.ct = ct,
                                verbose = TRUE)
prop_file = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
                  "prop_primary_tumor_unstranded_subset_CID3586.rds", 
                  sep="")
msg <- paste("Saving props to ", prop_file)
print(msg)
saveRDS(prop, file=prop_file)

#undebug(music_prop)