rm(list = ls())

library(ggplot2)
library(dplyr)


setwd("~/github/compbio_project/code")

 
sce_file = paste("../data/brca/GSE161529/processed/",
                 "tumor_sce.rds",
                 sep="")
tumor_patients <- readRDS(sce_file)
table(tumor_patients$cell_types)

sce_file = paste("../data/brca/GSE161529/processed/", 
                 "normal_sce.rds", 
                 sep="")
normal_patients <- readRDS(sce_file)
#normal_patients <- normal_patients[1:10, ]
table(normal_patients$cell_types)

df <- tumor_patients

# Extract expression values for marker genes
marker_gene_counts <- as.matrix(counts(df))
marker_gene_counts <- as.data.frame(marker_gene_counts)

# Join marker gene counts with colData
df <-merge(t(marker_gene_counts), colData(df), by.x = 0, by.y = 0)

#  Group by individual and cell type, calculate mean expression for each group
first_gene <- rownames(marker_gene_counts)[[1]]
last_gene <- rownames(marker_gene_counts)[dim(marker_gene_counts)[[2]]]
df_summary <- df %>%
  group_by(Individual, cell_types) %>%
  mutate(mean_expression = mean(c_across(first_gene:last_gene))) %>%
  distinct(Individual, cell_types, .keep_all = TRUE) %>%
  select(Individual, cell_types, mean_expression)

df_summary


# Create violin plot
temp_ggplot <- ggplot(df_summary, aes(x = cell_types, y = mean_expression, fill = cell_types)) +
  geom_violin() +  # Add violins
  geom_point() +  # Add points for individual data points
  stat_summary(fun = mean, geom = "point", size = 3, shape = 21, fill = "white") + # Add mean points
  theme_classic() + 
  theme(plot.background = element_rect(size = 1.5, color = "black"),
        panel.background = element_rect(size = 1.5, color = "black"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(face = 'bold'),
        axis.text.y = element_text(size = 14, face = 'bold'),        
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.75, "cm"),
        legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Cell Type", y = "Mean Expression") +
  ggtitle("GSE161529 Mean Expression by Cell Type") 

# plot
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 12,
      height = 5,
      useDingbats=F
    )
  }
temp_pdf_function(paste("../figures/GSE161529/", 
                        "mean_expression_cell_type.pdf",
                        sep=""))
print(temp_ggplot)
dev.off()
