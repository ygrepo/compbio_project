rm(list = ls())

library(SummarizedExperiment)
library(TCGAbiolinks)
library(EPIC)
setwd("~/github/compbio_project/code")

query.exp <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)
GDCdownload(
  query = query.exp,
  files.per.chunk = 100
)

brca.exp <- GDCprepare(
  query = query.exp, 
  directory = "GDCdata",
  save = FALSE, 
  #save.filename = "brcaExp.rda"
)

#brca.exp = load(file='brcaExp.rda')
exprs_data <- as.data.frame(assays(brca.exp)$unstranded)

out <- EPIC(bulk = exprs_data, reference = "TRef")

#write.table(exprs_data, file = "brcaExp.txt", quote = FALSE, sep='\t')

# get subtype information
infomation.subtype <- TCGAquery_subtype(tumor = "BRCA")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical") 

# Which samples are Primary Tumor
samples.primary.tumour <- brca.exp$barcode[brca.exp$shortLetterCode == "TP"]

# which samples are solid tissue normal
samples.solid.tissue.normal <- brca.exp$barcode[brca.exp$shortLetterCode == "NT"]