rm(list = ls())

library(TCGAbiolinks)
library(SummarizedExperiment)


setwd('/Users/yvesgreatti/github/compbio_project/data')

load("GBM.rda")

str(data)
colData(data)
rowData(data)
assayNames(data)

# Extract the gene expression data
gene_expr <- assay(data)

# View the dimensions of the gene expression matrix
dim(gene_expr)

# View the first few rows and columns of the gene expression matrix
head(gene_expr[, 1:5])
