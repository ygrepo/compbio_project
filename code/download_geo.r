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

library(GEOquery)
library(limma)
library(dplyr)
library(tidyr)
library(scater)
library(tensorflow)
library(cellassign)
library(Matrix)

setwd("~/github/compbio_project/code")

geonumber = "GSE143423"
# associateSymbolsToGenes <- TRUE
# verbose <- TRUE
# geneExpressionDF <- getGeneExpressionFromGEO(geonumber,  associateSymbolsToGenes, verbose)
# # Save the vector in rds format
# saveRDS(geneExpressionDF, "../data/GSE222450")
# 
gset <- getGEO(geonumber)

eset = gset[[1]]
eset = makeSummarizedExperimentFromExpressionSet(eset)
eset = as(eset, "SingleCellExperiment")


