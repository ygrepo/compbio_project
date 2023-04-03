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
                 "tumor_sce.rds", 
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

cell.types <- c("AV","BA","Fibroblast","HS","Immune","Vascular and lymphatic")
cell.types.colors <- c("#56B4E9", "#F7DC6F",  "#009E73",
             "#E74C3C", "#0072B2", "#E67E22")

temp_colour_pal_2 <- data.frame(celltype = cell.types, colour = cell.types.colors)
# Create the bar plot
temp_ggplot <- ggplot(data = count_by_individual, aes(x = Individual, y = prop, fill = cellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x=element_text(angle=45,
                                 size=17.5, hjust=1),
        strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text.x = element_text(size=25, face = "bold"),
        axis.text = element_text(size=17.5), 
        axis.title = element_text(size=25, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.75, 'cm'), #change legend key size
        legend.key.height = unit(0.75, 'cm'), #change legend key height
        legend.key.width = unit(0.75, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=20) #change legend text font size
  ) + 
  # scale_y_continuovus(position = "right") +
  xlab("Patient ID") +
  ylab("Cell Proportions")  +
  guides(fill=guide_legend(title="Cell Type")) +
  # coord_flip()
  scale_fill_manual(values = as.vector(temp_colour_pal_2$colour)) 
#  theme_classic()

temp_pdf_function(paste("../figures/GS161529/", 
                        "tumor_cell_subtype_proportions_patients.pdf",
                        sep=""))
print(temp_ggplot)
dev.off()