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
library(SingleCellExperiment)
library(DropletUtils)

setwd("~/github/compbio_project/code")

# Then, use the readMM() function to read the MTX file
path_to_data = paste("../data//brca/Wu_etal_2021_BRCA_bulk_sc/", 
                  "Wu_etal_2021_BRCA_scRNASeq/", 
                  sep="")

sce <- read10xCounts(path_to_data)

# Read in the metadata file as a data frame
metadata_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/Wu_etal_2021_BRCA_scRNASeq/", 
                     "metadata.csv", 
                     sep="")
metadata <- read.csv(metadata_file, row.names = 1)

# Assign the columns to the colData() of the SingleCellExperiment object
sce@colData$cellType <- metadata$celltype_major
ct <- unique(sce$cellType)
print(ct)

sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                      "sce_major_cell_type.rds", 
                      sep="")
saveRDS(sce, file=sce_file)

sce_subset <- sce[, 1:1000]
sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "sce_1000_major_cell_type.rds", 
                 sep="")
saveRDS(sce_subset, file=sce_file)