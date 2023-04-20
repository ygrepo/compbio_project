  rm(list = ls())
  
  # The easiest way to get ggplot2 is to install the whole tidyverse:
  #install.packages("tidyverse")
  
  
  # load
  library(SingleCellExperiment)
  library(ggplot2)
  #library(easyGgplot2)
  library(Matrix)
  
  setwd("~/github/compbio_project/code")
  
  prop_file = paste("../data/brca/tcga/processed/GSE176078/ct_minor/", 
                    "prop_primary_tumor_unstranded_subset_CID3586.rds", 
                    sep="")
  prop <- readRDS(file=prop_file)
  head(prop$Est.prop.weighted)
  head(prop$r.squared.full)
  head(prop$Var.prop)

# Convert the matrix to a data frame
weight_df <- as.data.frame(prop$Est.prop.weighted)
weight_df$Individual <- rownames(weight_df)
weight_df_long <- tidyr::gather(weight_df, "cellType", "prop", -Individual)

# Get unique individual names and create ID mapping
individuals <- unique(weight_df_long$Individual)
id_mapping <- paste0("ID-", match(individuals, individuals))

# Replace individual names with IDs in weight_df_long
weight_df_long$Individual <- id_mapping[match(weight_df_long$Individual, individuals)]

# Sample n individuals
n <- 100
sampled_individuals <- sample(unique(weight_df_long$Individual), n)

# Filter weight_df_long to only include the sampled individuals
sampled_weight_df_long <- weight_df_long[weight_df_long$Individual %in% sampled_individuals, ]

immune.cells <- c("B cells Memory",
                  "B cells Naive",
                  "T cells CD8+",
                  "T cells CD4+",
                  "NK cells",
                  "Cycling T-cells",
                  "NKT cells",
                  "Monocyte")
colours <- c("#56B4E9", "#F7DC6F",  "#009E73",
             "#E74C3C", "#0072B2", "#E67E22", 
             "#A569BD", "#0B5345")

# immune.cells <- c("B cells Memory",
#                   "B cells Naive",
#                   "T cells CD8+",
#                   "T cells CD4+",
#                   "NK cells",
#                   "Cycling T-cells",
#                   "NKT cells",
#                   "Macrophage",
#                   "Monocyte",
#                   "Other")
# colours <- c("#56B4E9", "#F7DC6F",  "#009E73",
#              "#E74C3C", "#0072B2", "#E67E22", "#641E16",
#              "#A569BD", "#0B5345","#17202A")
# colours <- c("#56B4E9", "#F0E442",  "#009E73",
#              "#E69F00", "#0072B2", "#D55E00", "#CE79A7",
#              "#9C79A7", "#299999","#299949")

temp_colour_pal_2 <- data.frame(celltype = immune.cells, colour = colours)
# Create the bar plot
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 12,
      height = 5,
      useDingbats=F
    )
  }

temp_ggplot <- ggplot(data = sampled_weight_df_long, aes(x = Individual, y = prop, fill = cellType)) +
  geom_bar(stat = "identity") +
  labs(y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x=element_text(angle=45,
                                 size=12, hjust=0.3),
        strip.background = element_rect(colour="black", fill="white",
                                        size=0.1, linetype="solid"),
        strip.text.x = element_text(size=12, face = "bold"),
        axis.text = element_text(size=12), 
        axis.title = element_text(size=14, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.75, 'cm'), #change legend key size
        legend.key.height = unit(0.75, 'cm'), #change legend key height
        legend.key.width = unit(0.75, 'cm'), #change legend key width
        legend.title = element_text(size=12, face="bold"), #change legend title font size
        legend.text = element_text(size=12), #change legend text font size
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_blank()
  ) + 
  ylab("Cell Proportions")  +
  guides(fill=guide_legend(title="Cell Type")) +
  scale_fill_manual(values = as.vector(temp_colour_pal_2$colour)) +
  theme(axis.text.x = element_blank())

temp_pdf_function(paste("../figures/Wu/With_No_Other/", 
                        "Tcga_Wu_cell_type_proportions_patients_CID3586.pdf",
                        sep=""))
print(temp_ggplot)
dev.off()


df <- as_tibble(weight_df_long)

# Remove rows with 0 values in the 'prop' column
df <- subset(df, prop != 0)
# Create violin plot
temp_ggplot <- ggplot(df, aes(x = cellType, y = -log(prop), fill = cellType)) +
  geom_violin() +
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
  labs(x = "Cell Type", y = "-Log(Cell Proportion)", fill = "Cell Type")

temp_pdf_function(paste("../figures/GSE176078/With_No_Other/", 
                        "Tcga_Wu_cell_type_proportions_patients_CID3586_violin.pdf",
                        sep=""))
print(temp_ggplot)
dev.off()

