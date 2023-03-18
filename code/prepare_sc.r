rm(list = ls())

#install.packages("geneExpressionFromGEO", repos='http://cran.us.r-project.org')
#BiocManager::install("scater")
# install.packages("tensorflow")
# library(tensorflow)
# install_tensorflow(extra_packages = "tensorflow-probability")
# BiocManager::install('cellassign')
# install.packages("Matrix")


# load
library("geneExpressionFromGEO")
library(SingleCellExperiment)

#library(GEOquery)
library(limma)
library(dplyr)
library(tidyr)
library(scater)
library(Matrix)

setwd("~/github/compbio_project/code")

path_to_data = "/Users/yvesgreatti/github/compbio_project/data/"

count_matrix <- read.csv(paste(path_to_data, 
                "GSE143423_tnbc_scRNAseq_gene_expression_counts.csv", 
                sep=","), row.names = 1)
                         
# Read in the metadata file as a data frame
metadata_list <- read.csv(paste(path_to_data, 
                              "GSE143423_tnbc_scRNAseq_metadata.csv", 
                              sep=""), row.names = 1)
names(metadata_list) <- c("cellID", "sampleID")

music_prop(bulk.mtx = exprs_data, sc.sce = sce, verbose = F)
