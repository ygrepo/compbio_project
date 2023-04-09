rm(list = ls())

library(SummarizedExperiment)
library(TCGAbiolinks)
setwd("~/github/compbio_project/code")

project = "TCGA-BRCA"
sample_type = c("Primary Tumor")
#sample_type = c("Solid Tissue Normal")
#sample_type = c("Primary Tumor","Solid Tissue Normal")

query.exp <- GDCquery(
  project = project, 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = sample_type
)
GDCdownload(
  query = query.exp,
  files.per.chunk = 100
)

tcga_rda_file = "tcga_brca.rds"
tcga_exp <- GDCprepare(
  query = query.exp, 
  directory = "GDCdata",
  save = FALSE, 
  save.filename = tcga_rda_file
)

# saveRDS(sub_tcga, file=tcga_rds_file)

#tcga_exp <- tcga_exp[,1:300]
exprs_data <- as.data.frame(assays(tcga_exp)$unstranded)
# tcga_file = paste("../data/brca/tcga/processed/",
#                   "primary_tumor_unstranded_300.rds", sep="")
tcga_file = paste("../data/brca/tcga/processed/", 
                   "normal_tissue_unstranded.rds", sep="")
saveRDS(exprs_data, file=tcga_file)


### More examples

# get subtype information
# information.subtype <- TCGAquery_subtype(tumor = "HNSC")

# get clinical data
# information.clinical <- GDCquery_clinic(project = "TCGA-HNSC",type = "clinical") 

# Which samples are Primary Tumor
# samples.primary.tumour <- tcga_exp$barcode[tcga_exp$shortLetterCode == "TP"]

# which samples are solid tissue normal
# samples.solid.tissue.normal <- tcga_exp$barcode[tcga_exp$shortLetterCode == "TP"]