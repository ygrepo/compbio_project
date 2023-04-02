rm(list = ls())

# The easiest way to get ggplot2 is to install the whole tidyverse:
#install.packages("tidyverse")


# load
library(SingleCellExperiment)
library(MuSiC)
library(DropletUtils)
library(ggplot2)
library(Matrix)

setwd("~/github/compbio_project/code")

prop_file = paste("../data/brca/tcga/processed/Wu/", 
                  "prop_primary_tumor_unstranded_subset_1_ct_minor.rds", 
                  sep="")
prop <- readRDS(file=prop_file)

sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                # "sce.rds", 
                "subset_1_ct_minor.rds", 
                 sep="")
sce <- readRDS(sce_file)


# Construct artificial bulk dataset. Use all 4 cell types: alpha, beta, gamma, delta
construct.full = bulk_construct(sce, 
                                clusters = 'cellType', 
                                samples = 'Barcode')
