rm(list = ls())

library(progeny)
library(ggplot2)
library(dplyr)
library(anndata)

setwd("~/github/compbio_project/code")


path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
             "prop_primary_tumor_unstranded_subset_CID3586.rds", 
             sep="")
cell_estimation<- readRDS(path)
head(cell_estimation)

props <- cell_estimation$Est.prop.weighted
props <- as.matrix(props)

tcga_file = paste("../data/brca/tcga/processed/", 
                  "converted_genes_primary_tumor_unstranded_exp_data.rds", 
                  sep="")
exprs <-readRDS(tcga_file)
#exprs <- exprs[,1:2]
msg <- paste("Tcga exp. dimensions", dim(exprs), sep=":")
print(msg)
# 
# tcga_genes <- toupper(rownames(exprs))
# rownames(exprs) <- tcga_genes
# 
# sce_genes <- toupper(rownames(weights))
# 
# 
# markers <- intersect(sce_genes, tcga_genes)
# msg <- paste("Markers", length(markers), sep=":")
# print(msg)
# msg <- paste("#Sce genes", length(sce_genes), sep=":")
# print(msg)
# msg <- paste("#Tcga genes", length(tcga_genes), sep=":")
# print(msg)

mat <- as.matrix(exprs)
row.names(mat) <- rownames(exprs)

gene_sce <- mat %*% props
dim(gene_sce)
#gene_sce[is.na(gene_sce) | is.infinite(gene_sce)] <- 0

# path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
#                   "genes_prop_cells_primary_tumor_unstranded_subset_CID3586.rds", 
#                   sep="")
# msg <- paste("Saving to ", path)
# print(msg)
# saveRDS(gene_sce, file=path)

# path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
#              "genes_prop_cells_primary_tumor_unstranded_subset_CID3586.h5ad", 
#              sep="")
# write_h5ad(AnnData(X=gene_sce), path)

path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
             "genes_prop_cells_primary_tumor_unstranded_subset_CID3586.csv", 
             sep="")
write.table(gene_sce, file=path, sep=",")
#write_h5ad(AnnData(X=gene_sce), path)

# Run PROGENy to predict the activity of transcription factors (TFs)
progeny_results <- progeny(gene_sce, organism="Human")



