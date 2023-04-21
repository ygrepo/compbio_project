rm(list = ls())

library(progeny)
library(ggplot2)
library(dplyr)
library(anndata)

setwd("~/github/compbio_project/code")


path = paste("../data/brca/tcga/processed/GSE161529/", 
             "prop_cell_type_normal_tissue_unstranded.rds", 
             sep="")
# path = paste("../data/brca/tcga/processed/GSE161529/", 
#              "prop_cell_type_tumor_tissue_unstranded.rds", 
#              sep="")

cell_estimation<- readRDS(path)
head(cell_estimation)

props <- cell_estimation$Est.prop.weighted
props <- as.matrix(props)

prop_patients <- unique(rownames(props))
print(length(prop_patients))


# define a function to extract the first three parts
extract_second_three <- function(input_string) {
  string_parts <- strsplit(input_string, "-")[[1]]
  selected_parts <- paste(string_parts[2:3], collapse = "-")
  return(selected_parts)
}


# apply the function to the list of strings
prop_patients <- lapply(prop_patients, extract_second_three)


tcga_file = paste("../data/brca/tcga/processed/", 
                  "normal_tissue_unstranded.rds", 
                  sep="")
# tcga_file = paste("../data/brca/tcga/processed/", 
#                   "converted_genes_primary_tumor_unstranded_exp_data.rds", 
#                   sep="")
exprs <-readRDS(tcga_file)
#exprs <- exprs[,1:2]
msg <- paste("Tcga exp. dimensions", dim(exprs), sep=":")
print(msg)

tcga_patients <- unique(toupper(colnames(exprs)))
print(length(tcga_patients))

tcga_patients <- lapply(tcga_patients, extract_second_three)
common_patients <- intersect(tcga_patients, prop_patients)
print(paste0("Found common patients:", length(common_patients)))

colnames(exprs) <-  tcga_patients
exprs <- exprs[, colnames(exprs) %in% common_patients]

rownames(props) <- prop_patients
props <- props[rownames(props) %in% common_patients, ]
print(dim(props))
unlisted <- unlist(rownames(props))
props <- props[!(rownames(props) %in% unlisted[duplicated(unlisted)]), ]
props <- as.matrix(props)


#sce_genes <- toupper(rownames(props))

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

path = paste("../data/brca/tcga/processed/GSE161529/", 
             "genes_prop_cells_normal.csv", 
             sep="")
write.table(gene_sce, file=path, sep=",")
#write_h5ad(AnnData(X=gene_sce), path)

# Run PROGENy to predict the activity of transcription factors (TFs)
progeny_results <- progeny(gene_sce, organism="Human")



