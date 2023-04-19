rm(list = ls())

#library(limma)
library(GEOquery)
library(oligo)
library(tidyverse)
library(purrr)
library(edgeR)
library(SingleCellExperiment)
library(Matrix)
library(knitr)

setwd("~/github/compbio_project/code")

gse161529 <- getGEO('GSE161529')
#getGEOSuppFiles("GSE161529")
gse161529 <- gse161529[[1]]
pd <- pData(gse161529)
print(unique(pd$`cancer type:ch1`))
pd <- pd[pd$`cancer type:ch1` != "Normal",]
# Define paths to the gene, barcode, and mtx files

# Define the directory path
#directory <- "./GSE161529/zipped/"

# Find all barcode files ending in "*barcodes.tsv.gz" in the directory
#barcode_files <- list.files(directory, pattern = "*barcodes.tsv.gz", full.names = TRUE)
#mtx_files <- list.files(directory, pattern = "*matrix.mtx.gz", full.names = TRUE)


path = "./GSE161529/SCP1731/metadata/metadataInfo.txt"
metadata <- read.table(path, header=TRUE, sep="\t", quote="", comment.char="")
metadata <- metadata[2:nrow(metadata),]
metadata <- metadata[!(metadata$donor_id %in% c("Human-WT-A",
                                              "Human-WT-B",
                                              "Human-WT-C",
                                              "Human-WT-D")), ]
print(unique(metadata$donor_id))
# Remove the numbers from the text column using gsub() and a regular expression
metadata$biosample_id <- gsub("^\\d+\\.\\s+", "", metadata$biosample_id)
metadata$cell_subtypes <- gsub("^\\d+\\.\\s+", "", metadata$cell_subtypes)
print(unique(metadata$biosample_id))
print(unique(metadata$cell_subtypes))

metadata_output <- metadata %>%
  dplyr::group_by(biosample_id, cell_subtypes) %>%
  dplyr::select(biosample_id, cell_subtypes) %>%
  dplyr::distinct() %>%
  dplyr::arrange(biosample_id, cell_subtypes) %>%
  kable(format = "html",
        col.names = c("Major Type", "Minor Type"),
        caption="GSE161529 Cell Type Hierarchy")

writeLines(metadata_output, "../figures/GSE161529/GSE161529_cell_hierarchy.html")


metadata$barcode <- gsub(".*_", "", gsub(".*-", "", metadata$NAME))
m_barcodes <- unique(metadata$barcode)
print(length(m_barcodes))

barcode_files <- str_split(pd$supplementary_file_1,"/") %>% map_chr(tail,1)
matrix_files <- str_split(pd$supplementary_file_2,"/") %>% map_chr(tail,1)
comb_files = list()
m_list = list()
comb_files$bc <- barcode_files
comb_files$m <- matrix_files
gene_files <- c("./GSE161529/GSE161529_features.tsv")
barcode_l = list()
gene_symbol_l = list()
biosample_id_l = list()
cell_types_l = list()
cell_subtypes_l = list()
donor_l = list()
for (i in (1:length(comb_files$bc))) {
  print(i)
  # Read the files using read10X
  barcode_file = paste0("./GSE161529/zipped/", comb_files$bc[[i]])
  matrix_file = paste0("./GSE161529/zipped/", comb_files$m[[i]])
  dge <- read10X(genes = gene_files, 
                  barcodes = barcode_file, 
                  mtx = matrix_file,
                  DGEList = FALSE)
  # Remove characters after "_" in the "barcode" column
  barcodes <- unique(gsub("-1", "", gsub("_.*", "", dge$samples$Barcode)))
  int_barcodes <- intersect(barcodes, m_barcodes)
  n_intersection <-length(int_barcodes)
  print(paste0("Intersection:", n_intersection))
  if (n_intersection == 0){
    print("Skipping...")
    next
  }
  barcode_l <- c(barcode_l,int_barcodes)
  colnames(dge$counts) <- barcodes
  gene_symbols <- toupper(dge$genes$Symbol)
  rownames(dge$counts) <- gene_symbols
  gene_symbol_l <- c(gene_symbol_l, gene_symbols)
  int_counts <- dge$counts[, colnames(dge$counts) %in% int_barcodes]
  print(paste0("Dim:", dim(int_counts)))
  m <- Matrix(int_counts)
  m_list <- c(m_list,m)
  s_metadata <- metadata[metadata$barcode %in% int_barcodes, ]
  biosample_id_l <- c(biosample_id_l, s_metadata$s_metadata)
  donor_l <- c(donor_l, s_metadata$donor_id)
  print(paste0("Length cell types:", length(s_metadata$biosample_id)))
  cell_types_l <- c(cell_types_l, s_metadata$biosample_id)
  print(paste0("Length cell subtypes:", length(s_metadata$cell_subtypes)))
  cell_subtypes_l <- c(cell_subtypes_l, s_metadata$cell_subtypes)
}

m <- do.call(cbind2, m_list)
# Create a SingleCellExperiment object with m
sce <- SingleCellExperiment(list(counts=m))
rowData(sce)$gene <-   do.call(cbind, gene_symbol_l)
sce$Individual <- do.call(cbind, donor_l)
sce$sampleID <- do.call(cbind, barcode_l)
sce$cell_types <- do.call(cbind, cell_types_l)
sce$cell_subtypes <- do.call(cbind, cell_subtypes_l)
sce$biosample_id <- do.call(cbind, biosample_id_l)
# Print the resulting SingleCellExperiment object
print(sce)


tcga_file = paste("../data/brca/tcga/processed/", 
                  "converted_genes_primary_tumor_unstranded_exp_data.rds", 
                  sep="")
exprs <-readRDS(tcga_file)

tcga_genes <- toupper(rownames(exprs))
sce_genes <- toupper(rowData(sce)$gene)

markers <- intersect(sce_genes, tcga_genes)
msg <- paste("Markers", length(markers), sep=":")
print(msg)
msg <- paste("#sce genes", length(sce_genes), sep=":")
print(msg)
msg <- paste("#tcga genes", length(tcga_genes), sep=":")
print(msg)

print(unique(sce$cell_types))

sce_file = paste("../data/brca/GSE161529/processed/", 
                 "tumor_sce.rds", 
                 sep="")
saveRDS(sce, file=sce_file)

# Filter the SingleCellExperiment object based on a metadata column
filter_values <- c("Immune")
immune_sce <- subset(sce, ,colData(sce)$cell_types == 'Immune')
sce_file = paste("../data/brca/GSE161529/processed/", 
                 "tumor_immune_sce.rds", 
                 sep="")
saveRDS(immune_sce, file=sce_file)


# Examples
# gse33146 <- getGEO('GSE33146')
# getGEOSuppFiles("GSE33146")
# gse33146 <- gse33146[[1]]
# pd <- pData(gse33146)
# gse33146_celdata <- read.celfiles(paste0('GSE33146/',pd$cel_file),
#                                   phenoData=phenoData(gse33146))
# setwd("GSE33146/")
# gse33146_celdata <- read.celfiles(dir(".", "CEL.gz"), phenoData=phenoData(gse33146))
# pData(gse33146_celdata)[,c("geo_accession","cell line:ch1","culture medium:ch1")]
# 
# gse66417 <- getGEO('GSE66417')
# getGEOSuppFiles("GSE66417")
# gse66417 <- gse66417[[1]]
# pd <- pData(gse66417)
# pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
# gse66417_celdata <- read.celfiles(paste0('GSE66417/',
#                                          pd$cel_file),
#                                   phenoData=phenoData(gse66417))
# pData(gse66417_celdata)[,c("geo_accession","cell type:ch1","treatment:ch1")]
# 
# gse65194 <- getGEO('GSE65194')
# getGEOSuppFiles("GSE65194")
# gse65194 <- gse65194[[1]]
# pd2 <- pData(gse65194)
# 
# 
# 
# gse1569698 <- getGEO('GSE156969')
# getGEOSuppFiles("GSE156969")
# gse1569698 <- gse1569698[[1]]
# pd3 <- pData(gse1569698)
# pd3['cel_file'] <- str_split(pd3$supplementary_file,"/") %>% map_chr(tail,1)
# gse164898_celdata <- read.celfiles(paste0('GSE164898/',
#                                          pd3$cel_file),
#                                   phenoData=phenoData(gse164898))


