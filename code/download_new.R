rm(list = ls())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("SummarizedExperiment")



library(TCGAbiolinks)
library(SummarizedExperiment)

setwd('/Users/yvesgreatti/github/compbio_project/data')

tcga_list = c('TCGA-LUAD')

for(i in tcga_list){
  ### expression ###
  query <- GDCquery(project = i,
                    data.category = 'Transcriptome Profiling',
                    data.type = 'Gene Expression Quantification',
                    workflow.type = 'STAR - Counts',
                    sample.type = 'Primary Tumor')
  if(is.null(query) != TRUE){
    GDCdownload(query)
    data <- GDCprepare(query)
    eset <- assay(data)
    write.csv(eset, file = paste0('/Users/yvesgreatti/github/compbio_project/data/GDCdata/',i,'/rna_seq_table.csv'))
  }
}
