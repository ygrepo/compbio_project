rm(list = ls())

# install.packages("tidyverse")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")

# library(tidyverse)
# library(data.table)
# library(edgeR)

# install.packages("remotes")
# remotes::install_github("GfellerLab/EPIC")

library("EPIC")

library(TCGAbiolinks)
library(SummarizedExperiment)

setwd('/Users/yvesgreatti/github/compbio_project/data')

query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "results", 
                  sample.type = c("Primary Tumor"),
                  legacy = TRUE)
# We will use only 10 samples to make the example faster
query$results[[1]] <-  query$results[[1]]                
#query$results[[1]] <-  query$results[[1]][1:10,]                  
GDCdownload(query)
gbm.exp <- GDCprepare(query, save = TRUE, summarizedExperiment = TRUE, save.filename = "GBMIllumina_HiSeq.rda")

load("GBMIllumina_HiSeq.rda")
gene_expression_matrix <- assay(gbm.exp)
gene_expression_df <- as.data.frame(gene_expression_matrix)
#gene_expression_df <- t(gene_expression_df)

#write.epic(gene_expression_df, file = "/Users/yvesgreatti/github/compbio_project/data/gbm.exp.txt", sep = "\t", quote = FALSE)

write.table(gene_expression_df, file = "/Users/yvesgreatti/github/compbio_project/data/tcga_gbm.exp.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

query <- GDCquery(project = "TCGA-GBM", data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM")
GDCdownload(query)
            
# Load gene expression data and sample information
exp_data <- GDCprepare(query, save=FALSE, directory = "path/to/downloaded/files")




# # Set the path to the top-level directory
# top_level_dir <- '/Users/yvesgreatti/github/compbio_project/data/GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification'
# 
# 
# # List all the sub-directories under the top-level directory
# sub_dirs <- list.dirs(path = top_level_dir, recursive = FALSE)
# 
# # Loop through each sub-directory and read in the gene count file
# gene_counts_list <- lapply(sub_dirs, function(sub_dir) {
#   
#   # List all the files in the sub-directory with .tsv extension
#   files <- list.files(path = sub_dir, pattern = "\\.tsv$", full.names = TRUE)
#   
#   # Read in the gene count file
#   gene_counts <- read.delim(files, header = TRUE, row.names = 1, sep = "\t")
#   
#   # Return the gene counts
#   return(gene_counts)
# })

