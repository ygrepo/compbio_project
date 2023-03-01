rm(list=ls())
# library(devtools)
# install_github('dviraran/xCell')
# install_github('GfellerLab/EPIC')
# install_github('ebecht/MCPcounter/Source')
# devtools::install_github('icbi-lab/immunedeconv')
library(immunedeconv)

#rna_file <- read.table("data/GSE81861_CRC_tumor_all_cells_FPKM.csv", header=TRUE, sep="\t")
# rna_file <- read.table("data/GSE81861/GSE81861_CRC_tumor_all_cells_FPKM.csv", header=FALSE, sep=",")
# eset <- new("ExpressionSet", exprs=as.matrix(rna_file))
data <- read.csv("data/GSE81861/GSE81861_CRC_tumor_all_cells_FPKM.csv", row.names=1)

data <- t(data)
exprs_matrix <- data.matrix(data)
eset <- new("ExpressionSet", exprs=exprs_matrix)
dm <- deconvolute(exprs_matrix, method="mcp_counter")

data(example_gene_expression_matrix)
timer_available_cancers
deconvolute(gene_expression, method=“timer”, indications = c (“BRCA”, “BRCA”, “LIHC”, “LIHC”))