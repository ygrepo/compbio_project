rm(list = ls())
setwd("~/github/OIO_Shaoshi/Deciphering_TME/Preprocess/")

path = "./data/brca/normal/scRNA/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct"
# Read in the .gct file
data <- read.table(path, header=TRUE, sep="\t", quote="", comment.char="", skip=2)

# Extract the gene expression matrix
gene_expression <- data[,-(1:2)]

# Extract the sample annotations
sample_annotations <- data.frame(data[,1:2])

# Remove unnecessary columns
gene_expression <- gene_expression[,3:ncol(gene_expression)]
