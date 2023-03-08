
rm(list = ls())

# install devtools if necessary
#install.packages('devtools')

# The following initializes usage of Bioc devel
#BiocManager::install(version = '3.16')

#BiocManager::install("TOAST")

# devtools::install_github('xuranw/MuSiC')

# load
library(MuSiC)
library(SingleCellExperiment)

setwd("~/github/compbio_project/code")

GSE50244.bulk.eset = readRDS('../data/GSE50244bulkeset.rds')
#GSE50244.bulk.eset

bulk.mtx = exprs(GSE50244.bulk.eset)
bulk.meta = exprs(GSE50244.bulk.eset)
EMTAB.sce = readRDS('../data/EMTABsce_healthy.rds')

#Estimate cell type proportions
Est.prop.GSE50244 = music_prop(bulk.mtx = bulk.mtx, sc.sce = EMTAB.sce, clusters = 'cellType',
                               samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma',
                                                                   'acinar', 'ductal'), verbose = F)

