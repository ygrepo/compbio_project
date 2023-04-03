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

dfs <- df %>%
  dplyr::select(B.cells.Memory, 
                SUBTYPE, 
                AJCC_PATHOLOGIC_TUMOR_STAGE,
                OS_MONTHS, 
                OS_STATUS) %>%
  dplyr::filter(B.cells.Memory >= 0)

dfs <- as_tibble(dfs)

print(dim(dfs))
print(length(dfs$OS_MONTHS))
print(length(dfs$OS_STATUS))


# Determine the row index to split the data into two groups
n <- nrow(dfs)
split_index <- floor(n / 2)


# Add B.cells.Memory.group based on the split index
dfs <- dfs %>%
  mutate(
    cell_type_group = ifelse(row_number() <= split_index, "low", "high")
  )

dfs <- dplyr::select(dfs, cell_type_group, SUBTYPE, 
                     AJCC_PATHOLOGIC_TUMOR_STAGE,
                     OS_MONTHS, 
                     OS_STATUS)
print(dim(dfs))

# Split the data into two groups based on B.cells.Memory.group
low_group <- filter(dfs, cell_type_group == "low")
high_group <- filter(dfs, cell_type_group == "high")

# Fit Cox proportional hazards models for each group
low_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ SUBTYPE,
                   data = low_group)
high_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ SUBTYPE, 
                    data = high_group)


# Print the model summaries
summary(low_model)
summary(high_model)

# Compute the survival curves for each group
low_surv <- survfit(low_model)
high_surv <- survfit(high_model)

# Combine low_surv and high_surv objects
surv_list <- list("Low Density" = low_surv, "High Density" = high_surv)

# Compute the p-value for the difference between the two curves
p_value <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ SUBTYPE, 
                          data = dfs)$pvalue


# Plot the combined survival curves
surv_plot <- survminer::ggsurvplot_combine(surv_list, 
                                           data = dfs,
                   xlab = "Time (months)",
                   ylab = "Survival probability",
                   ggtheme = theme(plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   axis.title = element_text(size = 14, face = "bold"),
                                   legend.title = element_text(size = 14, face = "bold"),
                                   legend.text = element_text(size = 12, face = "bold")))
surv_plot$plot <- surv_plot$plot + 
  ggtitle("Survival plot for B.cells memory") + 
  labs(subtitle = paste0("Log-rank test p = " ,round(p_value,2)))
surv_plot

 