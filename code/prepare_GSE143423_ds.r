rm(list = ls())

# load
library(SingleCellExperiment)
library(DropletUtils)

setwd("~/github/compbio_project/code")

sce_file = paste("../data/brca/GSE143423/", 
                      "GSE143423_tnbc_scRNAseq_gene_expression_counts.csv", 
                      sep="")
sce <- read.csv(sce_file, row.names = 1)

# Read in the metadata file as a data frame
metadata_file = paste("../data/brca/GSE143423/", 
                      "GSE143423_tnbc_scRNAseq_metadata.csv", 
                      sep="")
metadata <- read.csv(metadata_file, row.names = 1)
print(names(metadata))

# Assign the columns to the colData() of the SingleCellExperiment object
sce@colData$cellType <- metadata$celltype_minor
#sce@colData$cellType <- metadata$celltype_subset
sce@colData$celltypeMinor <- metadata$celltype_minor



# define individual IDs of interest
# individuals <- c("CID3586",
#                  "CID3921",
#                  "CID45171",
#                  "CID3838",
#                  "CID4066",
#                  "CID44041",
#                  "CID4465",
#                  "CID4495",
#                  "CID44971",
#                  "CID44991")

individuals <- c("CID3586")

                 # "CID3921",
                 # "CID3941")

# subset SingleCellExperiment object by individual IDs
subset_sce <- subset(sce, ,colData(sce)$Individual %in% individuals)

# define immune cells
immune.cells <- c("B cells Memory",
                 "B cells Naive",
                 "T cells CD8+",
                 "T cells CD4+",
                 "NK cells",
                 "Cycling T-cells",
                 "NKT cells",
                 "Macrophage",
                 "Monocyte"
)
#immune.cells <- c("B-cells", "T-cells")
#subset_sce <- sce
# cond = ((colData(subset_sce)$cellType %in% immune.cells) | 
#   (colData(subset_sce)$celltypeMinor == 'Macrophage') |
#   (colData(subset_sce)$celltypeMinor == 'Monocyte'))
cond = (colData(subset_sce)$cellType %in% immune.cells)
subset_sce <- subset(subset_sce, ,cond)
#t <- subset(sce, ,colData(sce)$celltypeMinor == 'Macrophage')

print(dim(sce))
print(dim(subset_sce))

ct <- unique(subset_sce$cellType)
print(ct)

sce_file = paste("../data/brca/Wu_etal_2021_BRCA_bulk_sc/processed/", 
                 "subset_1_ct_minor.rds", 
                 sep="")
saveRDS(subset_sce, file=sce_file)
