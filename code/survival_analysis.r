rm(list = ls())

# library(devtools)
# install_github("therneau/survival")
# BiocManager::install("cBioPortalData")
# install.packages("survminer")
# install.packages("webshot2")

#library(cBioPortalData)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
#library(lubridate)
#library(ggsurvfit)
#library(gtsummary)
#library(tidycmprsk)
#library(condsurv)
library(dplyr)
library(tidyr)
library(tibble)
#library(cluster)

setwd("~/github/compbio_project/code")

# The goal of this analysis is to identify potential prognostic factors for breast cancer 
# and to better understand the relationship between the immune system and cancer progression.\
# To perform survival analysis while controlling for confounding factors such as tumor stage, 
# one could include the confounding variable as a covariate in the Cox proportional hazards model. 
# This would allow us to examine the effect of the cell type of interest on survival 
# while adjusting for the effect of the confounding variable.

path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
                  "weights_cell_type_tumor_CID3586.csv", 
                  sep="")
cell_estimation = read.table(path,sep=',')
head(cell_estimation)

path = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
             "clinical_cell_type_tumor_CID3586.csv", 
             sep="")
survival = read.table(path,sep=',')
head(survival)
rownames(survival) <- survival$patientId

print(dim(cell_estimation))
print(dim(survival))

df <- merge(cell_estimation, survival,  by.x = 0, by.y = 0)
df$OS_STATUS <- as.numeric(substr(df$OS_STATUS, 1, 1))
df$DFS_STATUS <- as.numeric(substr(df$DFS_STATUS, 1, 1))
print(dim(df))
#df <- na.omit(df)
print(dim(df))

# Select the relevant variables and filter out negative values of the specified column
dfs <- df %>%
  dplyr::select(B.cells.Memory , 
                SUBTYPE, 
                AJCC_PATHOLOGIC_TUMOR_STAGE, 
                OS_MONTHS, 
                OS_STATUS) %>%
  dplyr::filter(B.cells.Memory  >= 0) %>%
  dplyr::arrange(desc(B.cells.Memory))


# Determine the row index to split the data into two groups
n <- nrow(dfs)
split_index <- floor(n / 2)
 
# Add cell_type_group based on the split index
dfs <- dfs %>%
   mutate(
     cell_type_group = ifelse(row_number() <= split_index, "high", "low"))
 
# dfs <- dfs %>%
#         dplyr::select(dfs, 
#                      cell_type_group, 
#                      OS_MONTHS,
#                      OS_STATUS)

# Fit Cox proportional hazards model for survival based on cell type group
surv_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cell_type_group, data = dfs)
summary(surv_model)

# Compute log-rank test p-value
logrank_pval <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cell_type_group, data = dfs)$pvalue

# Create Kaplan-Meier survival curves based on cell type group
km_surv <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cell_type_group, data = dfs)

# Plot the Kaplan-Meier curves with a legend and p-value
surv_plot <- survminer::ggsurvplot(km_surv, data = dfs, 
                        #pval = list(logrank_pval), 
                        #pval.method = TRUE,
                        legend.title = "Cell type group",
                        xlab = "Time (months)",
                        ylab = "Survival probability",
                        ggtheme = theme(plot.title = element_text(hjust = 0.5),
                                        plot.subtitle = element_text(hjust = 0.5),
                                        axis.title = element_text(size = 14, face = "bold"),
                                        legend.title = element_text(size = 14, face = "bold"),
                                        legend.text = element_text(size = 12, face = "bold")))

surv_plot

# Split the data into two groups based on the specified variable
low_group <- filter(data, cell_type_group == "low")
high_group <- filter(data, cell_type_group == "high")

# Fit Cox proportional hazards models for each group
low_model <- coxph(Surv(OS_MONTHS, OS_STATUS) , data = low_group)
high_model <- coxph(Surv(OS_MONTHS, OS_STATUS), data = high_group)

# Print the model summaries
cat("Low density group:\n")
print(summary(low_model))
cat("\nHigh density group:\n")
print(summary(high_model))


subtype_survival_analysis <- function(data, cell_type) {

  
  # Compute the survival curves for each group
  low_surv <- survfit(low_model)
  high_surv <- survfit(high_model)
  
  # Combine low_surv and high_surv objects
  surv_list <- list("Low Density" = low_surv, "High Density" = high_surv)
  
  # Compute the p-value for the difference between the two curves
  p_value <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = data)$pvalue
  
  # Plot the combined survival curves
  surv_plot <- survminer::ggsurvplot_combine(surv_list, 
                                             data = data,
                                             legend.title = cell_type,
                                             xlab = "Time (months)",
                                             ylab = "Survival probability",
                                             ggtheme = theme(plot.title = element_text(hjust = 0.5),
                                                             plot.subtitle = element_text(hjust = 0.5),
                                                             axis.title = element_text(size = 14, face = "bold"),
                                                             legend.title = element_text(size = 14, face = "bold"),
                                                             legend.text = element_text(size = 12, face = "bold")))
  surv_plot$plot <- surv_plot$plot + 
    ggtitle(paste0("Survival plot for ", cell_type)) + 
    labs(subtitle = paste0("Log-rank test p = " ,round(p_value,2)))
  
  return(surv_plot)
}

subtype_survival_analysis(dfs, "B.cells.Memory")
