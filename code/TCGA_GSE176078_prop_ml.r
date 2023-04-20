rm(list = ls())
  
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(tidyr)
library(dplyr)

setwd("~/github/compbio_project/code")
  
prop_file = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
                    "prop_primary_tumor_unstranded_subset_CID3586.rds", 
                    sep="")
prop <- readRDS(file=prop_file)
head(prop$Est.prop.weighted)
head(prop$r.squared.full)
head(prop$Var.prop)

weights <- prop$Est.prop.weighted
prop_patients <- unique(rownames(weights))
print(length(prop_patients))
#grep("^TCGA-A7-A13D", prop_patients, value = TRUE)


# define a function to extract the first three parts
extract_second_three <- function(input_string) {
  string_parts <- strsplit(input_string, "-")[[1]]
  selected_parts <- paste(string_parts[2:3], collapse = "-")
  return(selected_parts)
}


# apply the function to the list of strings
prop_patients <- lapply(prop_patients, extract_second_three)


estimate_file = paste("../data/brca/tcga/", 
                  "estimate.txt", 
                  sep="")
estimate_scores <- read.table(estimate_file,sep="\t", header = TRUE)
tcga_score_patients <- unique(estimate_scores$ID)
print(length(tcga_score_patients))

tcga_score_patients <- lapply(tcga_score_patients, extract_second_three)
common_patients <- intersect(tcga_score_patients, prop_patients)
print(paste0("Found common patients:", length(common_patients)))

rownames(weights) <- prop_patients
weights <- weights[rownames(weights) %in% common_patients, ]
print(dim(weights))
unlisted <- unlist(rownames(weights))
weights <- weights[!(rownames(weights) %in% unlisted[duplicated(unlisted)]), ]
weights <- as.data.frame(weights)


estimate_scores <- as.matrix(estimate_scores)
rownames(estimate_scores) <- tcga_score_patients
estimate_scores <- estimate_scores[rownames(estimate_scores) %in% common_patients, ]
print(dim(estimate_scores))
unlisted <- unlist(rownames(estimate_scores))
unlisted[duplicated(unlisted)]
estimate_scores <- estimate_scores[!(rownames(estimate_scores) %in% unlisted[duplicated(unlisted)]), ]
estimate_scores <- as.data.frame(estimate_scores)
estimate_scores <- select(estimate_scores, -ID)

# merge the data frames on row names
df_merged <- merge(weights, estimate_scores, by = "row.names")

# rename the row.names column to "id"
colnames(df_merged)[1] <- "ID"

df_merged$Immune_score <- as.numeric(as.character(df_merged$Immune_score))

df_merged$`B cells Memory`[df_merged$`B cells Memory` == 0] <- 1
df_merged$`B cells Memory` <- log(df_merged[,c("B cells Memory")])
df_merged$`B cells Naive`[df_merged$`B cells Naive` == 0] <- 1
df_merged$`B cells Naive` <- log(df_merged[,c("B cells Naive")])
df_merged$`T cells CD8+`[df_merged$`T cells CD8+` == 0] <- 1
df_merged$`T cells CD8+` <- log(df_merged[,c("T cells CD8+")])
df_merged$`T cells CD4+`[df_merged$`T cells CD4+` == 0] <- 1
df_merged$`T cells CD4+` <- log(df_merged[,c("T cells CD4+")])
df_merged$`NK cells`[df_merged$`NK cells` == 0] <- 1
df_merged$`NK cells` <- log(df_merged[,c("NK cells")])
df_merged$`NKT cells`[df_merged$`NKT cells` == 0] <- 1
df_merged$`NKT cells` <- log(df_merged[,c("NKT cells")])
df_merged$`Monocyte`[df_merged$`Monocyte` == 0] <- 1
df_merged$`Monocyte` <- log(df_merged[,c("Monocyte")])

path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
             "prop_estimate_scores_primary_tumor_unstranded_subset_CID3586.csv", 
             sep="")
write.table(df_merged, file=path, sep=",")


