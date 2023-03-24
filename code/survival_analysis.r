rm(list = ls())
# library(devtools)
# install_github("therneau/survival")
# BiocManager::install("cBioPortalData")
# install.packages("survminer")
# install.packages("webshot2")

library(cBioPortalData)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)
library(dplyr)


setwd("~/github/compbio_project/code")

tcga_file = paste("../data/brca/tcga/processed/", 
                  "clinical.rds", 
                  sep="")
clinical <- readRDS(file=tcga_file)


table(clinical$OS_STATUS)
table(clinical$SEX)
clinical$OS_MONTHS <- as.numeric(clinical$OS_MONTHS)
clinical$AGE <- as.numeric(clinical$AGE)

head(clinical[, c("OS_MONTHS", "OS_STATUS", "SEX")])

head(Surv(clinical$OS_MONTHS, as.numeric(substr(clinical$OS_STATUS, 1, 1))))

clinical$OS_STATUS <- as.numeric(substr(clinical$OS_STATUS, 1, 1))

survfit2(Surv(OS_MONTHS, OS_STATUS) ~ 1, data = clinical) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )


png(file = "../figures/clinic_brca_overall_surv.png")


survfit2(Surv(OS_MONTHS, OS_STATUS) ~ SEX, 
         data = clinical) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()
dev.off()

fit <-coxph(Surv(OS_MONTHS, OS_STATUS) ~ SEX, 
      data = clinical) 
#%>% 
t1 <- tbl_regression(fit, exp = TRUE) 

gt::gtsave(as_gt(t1), file = "../figures/clinic_brca_coxph_sex_tbl_reg.png")

# survfit2(Surv(OS_MONTHS, OS_STATUS) ~ AGE, 
#          data = clinical_brca) %>% 
#   ggsurvfit() +
#   labs(
#     x = "Days",
#     y = "Overall survival probability"
#   ) + 
#   add_confidence_interval() +
#   add_risktable()


# coxph(Surv(OS_MONTHS, OS_STATUS) ~ AGE, 
#              data = clinical_brca) %>% 
#   tbl_regression(exp = TRUE) 


fit <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ AJCC_PATHOLOGIC_TUMOR_STAGE, 
      data = clinical)
t1 <- tbl_regression(fit, exp = TRUE) 
gt::gtsave(as_gt(t1), file = "../figures/clinic_brca_coxph_stage_tbl_reg.png")


fit <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ RACE, 
      data = clinical)
t1 <- tbl_regression(fit, exp = TRUE) 
gt::gtsave(as_gt(t1), file = "../figures/clinic_brca_coxph_race_tbl_reg.png")


