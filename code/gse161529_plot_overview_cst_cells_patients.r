rm(list = ls())

# The easiest way to get ggplot2 is to install the whole tidyverse:
#install.packages("tidyverse")


# load
library(SingleCellExperiment)
library(DropletUtils)
library(ggplot2)
library(Matrix)

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
  cellType = sce$cell_subtypes
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

# View the resulting data frame
print(count_by_individual)

print(unique(count_by_individual$cellType))
print(length(unique(count_by_individual$cellType)))

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

cell.subtypes <-c("AP","I3 T cell","HSa","F3 Fibroblast","VL3 Pericyte",
                  "HSx","BL","HSb",                      
                  "BAb","VL2 Vascular endothelial",
                  "VL1 Lymphatic endothelial","I1 Myeloid cell",          
                  "BAa","Fx Fibroblast","I5 Plasma cell","F1 Fibroblast",            
                  "BAx","I2 NK cell")
cell.subtypes.colors <- rainbow(18)

temp_colour_pal_2 <- data.frame(celltype = cell.subtypes, colour = cell.subtypes.colors)
# Create the bar plot
temp_ggplot <- ggplot(data = count_by_individual, aes(x = Individual, y = prop, fill = cellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x=element_text(angle=45,
                                 size=12, hjust=1),
        strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text.x = element_text(size=18, face = "bold"),
        axis.text = element_text(size=12), 
        axis.title = element_text(size=16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=12) #change legend text font size
  ) + 
  # scale_y_continuovus(position = "right") +
  xlab("Patient ID") +
  ylab("Cell Proportions")  +
  guides(fill=guide_legend(title="Cell Type")) +
  # coord_flip()
  scale_fill_manual(values = as.vector(temp_colour_pal_2$colour)) 
#  theme_classic()

temp_pdf_function(paste("../figures/GS161529/", 
                        "normal_cell_subtype_proportions_patients.pdf",
                        sep=""))
print(temp_ggplot)
dev.off()
