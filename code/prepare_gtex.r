rm(list = ls())

# load
library(SingleCellExperiment)
library(DropletUtils)

setwd("~/github/compbio_project/code")

count_file = paste("../data/brca/gtex/", 
                      "counts.csv", 
                      sep="")
counts <- read.table(count_file,sep=",", header = TRUE)
expr <- counts[1:nrow(counts) , 2:ncol(counts)]

# Read in the metadata file as a data frame
annotation_file = paste("../data/brca/gtex/", 
                      "annotation.csv", 
                      sep="")
annotation <- read.table(annotation_file, header=TRUE, sep=',')
print(names(annotation))

sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(expr)),
  colData = annotation
)
rowData(sce) <- counts$X

immune.cells = unique(colData(sce)$CellType[grepl("^Immune", colData(sce)$CellType)])
cond = (colData(sce)$CellType %in% immune.cells)
subset_sce <- subset(sce, ,cond)

# Extract the description above the parenthesis for "Immune" values
subset_sce$CellType <- sub("^Immune \\((.*)\\).*", "\\1", 
                    subset_sce$CellType[grepl("^Immune", subset_sce$CellType)])
subset_sce$sampleID <- subset_sce$Barcode
subset_sce$Barcode <- NULL

print(dim(sce))
print(dim(subset_sce))

ct <- unique(subset_sce$CellType)
print(ct)

sce_genes <- toupper(rowData(subset_sce)$value)
rowData(subset_sce)$value <- sce_genes
rownames(subset_sce) <- sce_genes

sce_file = paste("../data/brca/gtex/processed", 
                        "subset_sce.rds", 
                        sep="")
saveRDS(subset_sce, file=sce_file)
