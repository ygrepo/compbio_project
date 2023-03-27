rm(list = ls())


setwd("~/github/compbio_project/code")

prop_file = paste("../data/brca/tcga/processed/Wu/", 
                  "prop_primary_tumor_unstranded_subset_1_ct_minor.rds", 
                  sep="")

prop <- readRDS(prop_file)