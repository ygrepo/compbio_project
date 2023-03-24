rm(list = ls())

library(TCGAbiolinks)

# Gene expression aligned against hg38
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  #workflow.type = "HTSeq - FPKM-UQ",
)
GDCdownload(query)
data <- GDCprepare(query)
data <- gbm.exp.harmonized
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

datatable(assay(data)[1:100,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

rowRanges(data)