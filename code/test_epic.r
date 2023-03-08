rm(list = ls())


# BiocManager::install("GEOquery")
# BiocManager::install("JinmiaoChenLab/cytofkit")
#devtools::install_github('dviraran/xCell')

library("EPIC")
library(tidyverse)
library(GEOquery)
library(xCell) 


# download and extract the gene expression data
data1 = getGEO("GSE182109")
datExpr1 = exprs (data1[[1]])

# read the gene expression matrix from a TXT file
dataExprs <- read.table("../data/GSE194213_Raw_gene_counts_matrix.txt", header = TRUE, row.names = 1, sep = "\t")

# check the dimensions of the expression matrix
dim(dataExprs)

dataExprs <- as.data.frame(dataExprs)
dataExprs <- dataExprs[,-1]

#estimate cell type proportions using xCell
results <- xCellAnalysis(dataExprs)

# exprs <- as.data.frame(exprs(gse[[1]]))
