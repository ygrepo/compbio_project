rm(list = ls())

# install devtools if necessary
#install.packages('devtools')

# The following initializes usage of Bioc devel
#BiocManager::install(version = '3.16')

#BiocManager::install("TOAST")

# devtools::install_github('xuranw/MuSiC')
#install.packages("remotes")
#remotes::install_github("MarioniLab/DropletUtils")

# load
#library(MuSiC)
library(SingleCellExperiment)
library(DropletUtils)

setwd("~/github/compbio_project/code")

# Then, use the readMM() function to read the MTX file
path_to_data = "~/github/compbio_project/data/Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq"
sce <- read10xCounts(path_to_data)

# Set the column names in the metadata to SampleId
#colnames(sce@colData) <- c("sampleID", "Barcode")
#colnames(sce) <- colData(sce)$Barcode

# Read in the metadata file as a data frame
metadata_file = paste("../data/Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq/", "metadata.csv", sep="")
metadata <- read.csv(metadata_file, row.names = 1)

# Assign the columns to the colData() of the SingleCellExperiment object
sce@colData$cellType <- metadata$celltype_major

# sce_file = paste("../data/", "wu_etal_sce.rds", sep="")
# #saveRDS(sce, file=sce_file)
# sce <- readRDS(sce_file)
# 
# ct <- unique(sce$cellType)
# sce_subset <- sce[, 1:1000]
#assay(sce_subset, "counts") <- as.matrix(counts(sce_subset))
# sce_file = paste("../data/", "small_wu_etal_sce.rds", sep="")
# saveRDS(sce_subset, file=sce_file)

