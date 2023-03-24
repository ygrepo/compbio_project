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
library(MuSiC)
library(SingleCellExperiment)

setwd("~/github/compbio_project/code")
# source("music/utils.R")
# source("music/construct.R")
# source("music/analysis.R")


bulk_rna_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/",
                      "GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt",
                      sep="")
bulk_rna <- read.csv(bulk_rna_file, sep="\t")

# Select the expression data from the bulk RNA-seq data frame
bulk_expr <- as.matrix(bulk_rna[, -1])
bulk_expr <- bulk_expr [,1:2]
rownames(bulk_expr) <- bulk_rna$Genes

sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "sce_1000_major_cell_type.rds", 
                 sep="")
sce <- readRDS(sce_file)
ct <- unique(sce$cellType)


#debug(music_prop)
# Estimate cell type proportions
prop = music_prop(bulk.mtx = bulk_expr, sc.sce = sce, 
           clusters = 'cellType', 
           samples = 'Barcode',
           select.ct = ct,
           verbose = T)
prop_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "prop_2_1000.rds", 
                 sep="")
saveRDS(prop, file=prop_file)

# undebug(music_prop)

