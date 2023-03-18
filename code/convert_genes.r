rm(list = ls())

# install devtools if necessary
#install.packages('devtools')

# The following initializes usage of Bioc devel
#BiocManager::install(version = '3.16')

# BiocManager::install("EnsDb.Hsapiens.v79")

# devtools::install_github('xuranw/MuSiC')
#install.packages("remotes")
#remotes::install_github("MarioniLab/DropletUtils")
# BiocManager::install("AnnotationDbi", force = T)

#install.packages("dplyr")   

# load
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
#library(MuSiC)
library(SingleCellExperiment)
library(dplyr)

setwd("~/github/compbio_project/code")

source("music/utils.R")
source("music/construct.R")
source("music/analysis.R")


# tcga_file = paste("../data/", "small_brca_unstranded.rda", sep="")
# exprs_data <-readRDS(tcga_file)
# exprs_data <- exprs_data[,1:25]
# tcga_file = paste("../data/", "small_brca_unstranded.rda", sep="")
# saveRDS(exprs_data, file=tcga_file)

# Extract gene IDs from row names of exprs_data
# ensembl.genes <- sub('\\.[0-9]*$', '', rownames(exprs_data))
# 
# # Set row names of exprs_data to gene IDs
# rownames(exprs_data) <- ensembl.genes
# 
# # Select gene IDs and symbols from EnsDb object
# gene_ids <- ensembldb::select(EnsDb.Hsapiens.v79, 
#                               keys = ensembl.genes, keytype = "GENEID", 
#                               columns = c("SYMBOL", "GENEID"))
# bind_rows(gene_ids)
# 
# # Filter exprs_data by matching row names
# filtered_data <- exprs_data[row.names(exprs_data) %in% gene_ids$GENEID, ]
# #filtered_data <- head(filtered_data)
# filtered_data <- tibble::rownames_to_column(filtered_data, "row_names") 
# filtered_data <- rename(filtered_data, GENEID = row_names)
# filtered_data <- merge(filtered_data,gene_ids,by = "GENEID")
# filtered_data <-distinct(filtered_data, SYMBOL, .keep_all=TRUE)
# rownames(filtered_data) <- filtered_data[,ncol(filtered_data)]
# filtered_data <- filtered_data %>% select(-one_of(c('GENEID')))
# filtered_data <- filtered_data[,-ncol(filtered_data)]

tcga_file = paste("../data/", "filtered_brca_unstranded.rda", sep="")
#saveRDS(filtered_data, file=tcga_file)
filtered_data <-readRDS(tcga_file)


exprs <- as.matrix(filtered_data[, -1])
exprs <- exprs[,1:2]
rownames(exprs) <- rownames(filtered_data)


# 2. Convert from gene.symbol to ensembl.gene
#geneSymbols <-  c('IDS')
#geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, 
# keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

# sce_file = paste("../data/", "wu_etal_sce.rds", sep="")
# #saveRDS(sce, file=sce_file)
# sce <- readRDS(sce_file)


sce_file = paste("../data/", "small_wu_etal_sce.rds", sep="")
sce_subset <- readRDS(file=sce_file)
ct <- unique(sce_subset$cellType)



#debug(music_prop)
# Estimate cell type proportions
# debug(music_prop)
Est.prop.TCGA_prop = music_prop(bulk.mtx = exprs, sc.sce = sce_subset, 
                                clusters = 'cellType', 
                                samples = 'Barcode',
                                select.ct = ct,
                                verbose = TRUE)
prop_cells_file = paste("../data/", "031822_25_1000_prop.rds", sep="")
saveRDS(prop_cells, file=Est.prop.TCGA_prop)

#undebug(music_prop)