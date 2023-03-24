rm(list = ls())
# library(devtools)
# install_github("therneau/survival")
# BiocManager::install("cBioPortalData")
# install.packages("survminer")
# install.packages("webshot2")

# library(SummarizedExperiment)
#library(TCGAbiolinks)
library(cBioPortalData)
library(dplyr)


setwd("~/github/compbio_project/code")

# tcga_file = paste("../data/", "brca_300.rds", sep="")
# tcga <- readRDS(tcga_file)
# 

# tcga_file = paste("../data/brca_tcga_pan_can_atlas_2018/", "data_clinical_patient.txt", sep="")
# clinical_data <- read.csv(tcga_file, sep="\t")

cbio <- cBioPortal()
# # studies <- getStudies(cbio, buildReport = TRUE)
# cache_dir = "../data/cache"
# setCache(
#   directory = cache_dir,
#   verbose = TRUE,
#   ask = interactive()
# )


studyId = "brca_tcga_pan_can_atlas_2018" 
clinical <- clinicalData(cbio, studyId)
#clinical_file = paste("../data/", "brca_clinical_300.rds", sep="")

table(clinical$OS_STATUS)
table(clinical$SEX)

tcga_file = paste("../data/brca/tcga/processed/", 
                  "clinical.rds", 
                  sep="")
saveRDS(clinical, file=tcga_file)

