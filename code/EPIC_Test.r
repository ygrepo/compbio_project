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

setwd('/Users/yg/code/github/compbio_project/data')


# query <- GDCquery(project = "TCGA-LUAD", 
#                   legacy = TRUE,
#                   data.category = "Transcriptome Profiling", 
#                   data.type = "Gene Expression Quantification",
#                   workflow.type = "HTSeq - FPKM")
query <- GDCquery(project = c("TCGA-LUAD"),
                  data.category = "Sequencing Reads",  
                  sample.type = "Primary Tumor")

GDCdownload(query)
            
# Load gene expression data and sample information
exp_data <- GDCprepare(query, save=FALSE, directory = "data/TCGA-LUAD")




