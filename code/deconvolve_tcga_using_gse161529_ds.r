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

# tcga_file = paste("../data/brca/tcga/processed/", 
#                   "converted_genes_normal_tissue_unstranded_exp_data.rds", 
#                   sep="")
tcga_file = paste("../data/brca/tcga/processed/",
                  "converted_genes_primary_tumor_unstranded_exp_data.rds",
                  sep="")
exprs <-readRDS(tcga_file)
#exprs <- exprs[,1:2]
msg <- paste("Tcga exp. dimensions", dim(exprs), sep=":")
print(msg)

sce_file = paste("../data/brca/GSE161529/processed/", 
                 "tumor_immune_sce.rds", 
                 # "normal_sce.rds", 
                 sep="")
msg <- paste("Loading",
             sce_file,
             sep=" ")
print(msg)
sce <- readRDS(sce_file)
# sce_file = paste("../data/brca/GSE161529/processed/", 
#                  "tumor_immune_sce.rds", 
#                  sep="")

# immune_sce <- subset(sce, ,colData(sce)$cell_types == 'Immune')
# sce_file = paste("../data/brca/GSE161529/processed/", 
#                  "normal_immune_sce.rds", 
#                  sep="")
# saveRDS(immune_sce, file=sce_file)

msg <- paste("Sce dimensions", dim(sce), sep=":")
print(msg)

print(unique(sce$cell_types))
print(unique(sce$cell_subtypes))

#sce <- sce[,1:1000]
ct <- unique(sce$cell_subtypes)
#ct <- unique(sce$cell_subtypes)
msg <- paste("Cell types", ct, sep=":")
print(msg)
cst <- unique(sce$cell_subtypes)
msg <- paste("Cell subtypes", cst, sep=":")
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
                                clusters = 'cell_subtypes', 
                                samples = 'sampleID',
                                select.ct = ct,
                                verbose = TRUE)
print(head(prop$Est.prop.weighted))

prop_file = paste("../data/brca/tcga/processed/GSE161529/",
                  "prop_cell_subtype_immune_tumor_tissue.rds",
#                  "prop_cell_subtype_tumor_tissue_unstranded.rds",                  
                  sep="")
msg <- paste("Saving props to ", prop_file)
print(msg)
saveRDS(prop, file=prop_file)

#undebug(music_prop)