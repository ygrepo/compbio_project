rm(list = ls())

# The easiest way to get ggplot2 is to install the whole tidyverse:
#install.packages("tidyverse")


# load
library(SingleCellExperiment)
library(DropletUtils)
library(ggplot2)
library(Matrix)
library(tidyr)
library(dplyr)

setwd("~/github/compbio_project/code")

sce_file = paste("../data/brca/GSE161529/processed/", 
                 "normal_sce.rds", 
                 sep="")
sce <- readRDS(sce_file)

# Check the result
table(sce$cell_types)
table(sce$cell_subtypes)

# Compute the column sums of the sparse counts assay
col_sums <- colSums(assay(sce))

sce2 <- SingleCellExperiment(
  assays = list(counts = as.matrix(col_sums)),
)

# Create the rowData of sce2 using Individual 
# and cellType information from colData
rowData(sce2) <- DataFrame(
  Individual = sce$Individual,
  cellType = sce$cell_types
)

# Aggregate the counts by Individual and cellType
count_by_individual <- aggregate(x = assay(sce2), 
                                 by = list(Individual = rowData(sce2)$Individual, 
                                           cellType = rowData(sce2)$cellType), 
                                 FUN = sum)

# Compute the total count for each Individual
total_counts <- aggregate(x = assay(sce2), 
                          by = list(Individual = rowData(sce2)$Individual), 
                          FUN = sum)

# Merge the count and total count data frames
count_by_individual <- merge(count_by_individual, total_counts, by = "Individual")

# Compute the cell type proportion for each Individual
count_by_individual$prop <- count_by_individual$V1.x / count_by_individual$V1.y

count_by_individual <- select(count_by_individual, -V1.x)
count_by_individual <- select(count_by_individual, -V1.y)

# View the resulting data frame
print(count_by_individual)

print(unique(count_by_individual$cellType))
print(length(unique(count_by_individual$cellType)))

normal_df <- count_by_individual %>%
  pivot_wider(names_from = cellType, values_from = prop)
# group by Individual
normal_df <- normal_df %>%
  group_by(Individual) 
normal_df <- normal_df %>% ungroup() %>% select(colnames(prop$Est.prop.weighted))


prop_file = paste("../data/brca/tcga/processed/GSE161529/",
                  "prop_cell_type_tumor_tissue_unstranded.rds",
                  sep="")
# prop_file = paste("../data/brca/tcga/processed/GSE161529/",
#                   "prop_cell_type_normal_tissue_unstranded.rds",
#                   sep="")
msg <- paste("Reading props from ", prop_file)
print(msg)
prop <- readRDS(prop_file)


prop_est <- as.data.frame(prop$Est.prop.weighted)%>% select(colnames(prop$Est.prop.weighted))

head(prop_est)

print(wilcox.test(as.matrix(normal_df), as.matrix(prop_est)))


