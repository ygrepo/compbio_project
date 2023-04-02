rm(list = ls())

library(cBioPortalData)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)
library(dplyr)
library(tibble)

setwd("~/github/compbio_project/code")

tcga_file = paste("../data/brca/tcga/processed/", 
                  "clinical.rds", 
                  sep="")
clinical <- readRDS(file=tcga_file)


table(clinical$OS_MONTHS)
table(clinical$OS_STATUS)
table(clinical$DFS_MONTHS)
table(clinical$DFS_MONTHS)

tcga_clin_patients <- unique(clinical$patientId)
print(length(tcga_clin_patients))
#grep("^TCGA-A7-A13D", tcga_clin_patients, value = TRUE)

prop_file = paste("../data/brca/tcga/processed/GSE161529/", 
                  "prop_cell_type_tumor_tissue_unstranded.rds", 
                  sep="")
prop <- readRDS(file=prop_file)
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
tcga_clin_patients <- lapply(tcga_clin_patients, extract_second_three)
common_patients <- intersect(tcga_clin_patients, prop_patients)
print(paste0("Found common patients:", length(common_patients)))
rownames(weights) <- prop_patients
weights <- weights[rownames(weights) %in% common_patients, ]
print(dim(weights))
unlisted <- unlist(rownames(weights))
unlisted[duplicated(unlisted)]
weights <- weights[!(rownames(weights) %in% unlisted[duplicated(unlisted)]), ]

prop_file = paste("../data/brca/tcga/processed/GSE161529/", 
                  "weights_cell_type_tumor.csv", 
                  sep="")
write.table(weights,
            file=prop_file,
            quote=FALSE,
            sep=",",
            )
clinical$patientId <- tcga_clin_patients
clinical_tmp <- clinical[clinical$patientId %in% common_patients, ]
clinical_tmp <- clinical_tmp %>% select(patientId, OS_MONTHS, OS_STATUS, DFS_MONTHS, DFS_STATUS)
clinical_tmp$patientId <- vapply(clinical_tmp$patientId, paste, collapse = ", ", character(1L))
clin_file = paste("../data/brca/tcga/processed/GSE161529/", 
                  "clinical_cell_type_tumor.csv", 
                  sep="")
write.table(clinical_tmp,
            file=clin_file,
            quote=FALSE,
            sep=",",
)
