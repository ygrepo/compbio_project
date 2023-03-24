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
path_to_data = "../data/brca/Bhupinder_etal_2021_BRCA_bulk_sc/sc"
sce <- read10xCounts(path_to_data)
