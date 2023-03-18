rm(list = ls())

library(SummarizedExperiment)
library(TCGAbiolinks)
library(MuSiC)
library(SingleCellExperiment)

setwd("~/github/compbio_project/code")


path = "/Users/yvesgreatti/github/compbio_project/data/"

tcga_file = paste("../data/", "brca_unstranded.rda", sep="")
#save(exprs_data, file=tcga_rda_file)
exprs_data <- readRDS(file=tcga_file)
dim(exprs_data)
exprs_data <- exprs_data[,c(1,2)]

path_to_data = "/Users/yvesgreatti/github/compbio_project/data/"
count_matrix <- read.csv(paste(path_to_data, 
                               "GSE143423_tnbc_scRNAseq_gene_expression_counts.csv", 
                               sep=""), row.names = 1)

music_prop(bulk.mtx = exprs_data, sc.sce = sce, verbose = F)