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
# path_to_data = paste("../data//brca/Wu_etal_2021_BRCA_bulk_sc/",
#                   "Wu_etal_2021_BRCA_scRNASeq/",
#                   sep="")
# 
# sce <- read10xCounts(path_to_data)


# Read in the metadata file as a data frame
# metadata_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/Wu_etal_2021_BRCA_scRNASeq/", 
#                       "metadata.csv", 
#                       sep="")
# metadata <- read.csv(metadata_file, row.names = 1)
# print(unique(metadata$subtype))
# print(unique(metadata$celltype_subset))
# print(unique(metadata$celltype_minor))
# print(unique(metadata$celltype_major))
# 
# # Assign the columns to the colData() of the SingleCellExperiment object
# sce@colData$cellType <- metadata$celltype_minor
# #sce@colData$cellType <- metadata$celltype_subset
# sce@colData$celltypeMinor <- metadata$celltype_minor

# Define the immune cell types
# immune.cells <- c("B cells Memory",
#                   "B cells Naive",
#                   "T cells CD8+",
#                   "T cells CD4+",
#                   "NK cells",
#                   "Cycling T-cells",
#                   "NKT cells",
#                   "Macrophage",
#                   "Monocyte")

# Replace non-immune cell types with "other"
# colData(sce)$cellType <- ifelse(colData(sce)$cellType %in% immune.cells,
#                                 colData(sce)$cellType,
#                                 "Other")
# 
# # Check the result
# table(colData(sce)$cellType)
# 
# colData(sce)$Individual <- substring(colData(sce)$Barcode,
#                                      1,
#                                      regexpr("_", colData(sce)$Barcode)-1)
# 
sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/",
                 "sce.rds",
                 sep="")
#saveRDS(sce, file=sce_file)
sce <- readRDS(sce_file)


# define individual IDs of interest
individuals <- c("CID45171")
                 # "CID3586"
                 # "CID3963"
                 # "CID45171"

# subset SingleCellExperiment object by individual IDs
subset_sce <- subset(sce, ,colData(sce)$Individual %in% individuals)

#immune.cells <- c("B-cells", "T-cells")
#subset_sce <- sce
# cond = ((colData(subset_sce)$cellType %in% immune.cells) | 
#   (colData(subset_sce)$celltypeMinor == 'Macrophage') |
#   (colData(subset_sce)$celltypeMinor == 'Monocyte'))
#cond = (colData(subset_sce)$cellType %in% immune.cells)
#subset_sce <- subset(subset_sce, ,cond)

print(dim(sce))
print(dim(subset_sce))

ct <- unique(subset_sce$cellType)
print(ct)

sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "subset_CID45171.rds", 
                 sep="")
saveRDS(subset_sce, file=sce_file)
