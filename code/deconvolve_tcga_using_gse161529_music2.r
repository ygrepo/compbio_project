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
case_exprs <-readRDS(tcga_file)
#exprs <- exprs[,1:2]
msg <- paste("Tcga case exp. dimensions", dim(case_exprs), sep=":")
print(msg)

tcga_file = paste("../data/brca/tcga/processed/", 
                  "converted_genes_normal_tissue_unstranded_exp_data.rds", 
                  sep="")
control_exprs <-readRDS(tcga_file)
#exprs <- exprs[,1:2]
msg <- paste("Tcga control exp. dimensions", dim(control_exprs), sep=":")
print(msg)

sce_file = paste("../data/brca/GSE161529/processed/", 
                 "normal_sce.rds", 
                 sep="")
msg <- paste("Loading",
             sce_file,
             sep=" ")
print(msg)
sce <- readRDS(sce_file)
msg <- paste("Sce dimensions", dim(sce), sep=":")
print(msg)

ct <- unique(sce$cell_types)
msg <- paste("Cell types", ct, sep=":")
print(msg)

tcga_genes <- toupper(rownames(case_exprs))
rownames(case_exprs) <- tcga_genes

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
# music2 deconvolution
set.seed(1234)
est = music2_prop_t_statistics(bulk.control.mtx = control_exprs, 
                               bulk.case.mtx = case_exprs, 
                               sc.sce = sce, 
                               clusters = 'cell_types', 
                               samples = 'sampleID', 
                               select.ct = ct, 
                               n_resample=20, 
                               sample_prop=0.5,
                               cutoff_c=0.05,
                               cutoff_r=0.01)

est.prop = est$Est.prop
# prop = music_prop(bulk.mtx = case_exprs,
#                                 markers = markers,
#                                 sc.sce = sce, 
#                                 clusters = 'cell_types', 
#                                 samples = 'sampleID',
#                                 select.ct = ct,
#                                 verbose = TRUE)
prop_file = paste("../data/brca/tcga/processed/GSE161529/",
                  "prop_cell_type_music2_unstranded.rds",
                 sep="")
msg <- paste("Saving props to ", prop_file)
print(msg)
saveRDS(prop, file=prop_file)

#undebug(music_prop)