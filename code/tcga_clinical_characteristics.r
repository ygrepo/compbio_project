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
                  "clinical.rds", 
                  sep="")
tumor_clinical <- readRDS(file=tcga_file)

tumor_clinical$AGE <- as.numeric(as.character(tumor_clinical$AGE))  # Convert one variable to numeric

colnames(tumor_clinical)

summary(tumor_clinical$AGE)
sd(tumor_clinical$AGE)

Hmisc::describe(tumor_clinical$RACE)
Hmisc::describe(tumor_clinical$ETHNICITY)
Hmisc::describe(tumor_clinical$CANCER_TYPE)
Hmisc::describe(tumor_clinical$SUBTYPE)
Hmisc::describe(tumor_clinical$TUMOR_TYPE)
Hmisc::describe(tumor_clinical$AJCC_PATHOLOGIC_TUMOR_STAGE)
Hmisc::describe(tumor_clinical$RADIATION_THERAPY)
Hmisc::describe(tumor_clinical$RAGNUM_HYPOXIA_SCORE)
Hmisc::describe(tumor_clinical$WINTER_HYPOXIA_SCORE)

tcga_file = paste("../data/brca/tcga/processed/", 
                  "normal_tissue.rds", sep="")
normal_tcga <- readRDS(tcga_file)
normal_clinical <- colData(normal_tcga)
sapply(normal_clinical, class)
unique(normal_clinical$tumor_descriptor)
unique(normal_clinical$primary_site)
unique(normal_clinical$vital_status)
unique(normal_clinical$days_to_diagnosis)

Hmisc::describe(normal_clinical$age_at_index)

age_at_index <- normal_clinical %>% 
  as.data.frame() %>% 
  select(age_at_index)
summary(age_at_index)
sd(normal_clinical$age_at_index)

Hmisc::describe(normal_clinical$race)
Hmisc::describe(normal_clinical$ethnicity)
