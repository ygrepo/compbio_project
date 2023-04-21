rm(list = ls())
  
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(tidyr)
library(dplyr)

setwd("~/github/compbio_project/code")

rm(list = ls())

# library(cBioPortalData)
# library(survival)
# library(lubridate)
# library(ggsurvfit)
# library(gtsummary)
# library(tidycmprsk)
# library(condsurv)
library(dplyr)
library(tibble)
library(Hmisc)
library(knitr)

setwd("~/github/compbio_project/code")

tcga_file = paste("../data/brca/tcga/processed/", 
                  "normal_tissue.rds", sep="")
normal_tcga <- readRDS(tcga_file)
normal_patients <- unique(normal_tcga$patient)

# define a function to extract the first three parts
extract_second_three <- function(input_string) {
  string_parts <- strsplit(input_string, "-")[[1]]
  selected_parts <- paste(string_parts[2:3], collapse = "-")
  return(selected_parts)
}


# apply the function to the list of strings
normal_patients <- lapply(normal_patients, extract_second_three)
print(length(print(length(normal_patients))))


estimate_file = paste("../data/brca/tcga/", 
                  "estimate.txt", 
                  sep="")
estimate_scores <- read.table(estimate_file,sep="\t", header = TRUE)
tcga_score_patients <- unique(estimate_scores$ID)
print(length(tcga_score_patients))

tcga_score_patients <- lapply(tcga_score_patients, extract_second_three)
common_patients <- intersect(tcga_score_patients, normal_patients)
print(paste0("Found common patients:", length(common_patients)))


estimate_scores <- as.matrix(estimate_scores)
rownames(estimate_scores) <- tcga_score_patients
estimate_scores <- estimate_scores[rownames(estimate_scores) %in% common_patients, ]
print(dim(estimate_scores))
unlisted <- unlist(rownames(estimate_scores))
unlisted[duplicated(unlisted)]
estimate_scores <- estimate_scores[!(rownames(estimate_scores) %in% unlisted[duplicated(unlisted)]), ]
estimate_scores <- as.data.frame(estimate_scores)
estimate_scores <- select(estimate_scores, -ID)

path = paste("../data/brca/tcga/processed/", 
             "estimate_scores_normal_tissue.csv", 
             sep="")
write.table(estimate_scores, file=path, sep=",")


