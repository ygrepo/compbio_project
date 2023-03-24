rm(list = ls())

#BiocManager::install("SCAN.UPC")

library(SCAN.UPC)
library(GEOquery)

getGEOSuppFiles("GSE14520", makeDirectory=FALSE, baseDir=getwd())
tmpDir = "./data/GSE14520_RAW/"
celFilePath = file.path(tmpDir, "GSE14520.CEL.gz")
celFilePath = file.path(tmpDir, "GSM555237.CEL.gz")
celFilePath = file.path(tmpDir, "GSM363205.CEL.gz")

#normalized = SCAN("./data/GSE14520_RAW/GSM*")
normalized = SCAN(celFilePath)


getGEOSuppFiles("GSE115469", makeDirectory=FALSE, baseDir=getwd())
tmpDir = "./data/"
celFilePath = file.path(tmpDir, "GSE115469_Data.csv")
# Extract the expression data
data <- read.table(celFilePath, sep=",")
  
# Write the expression data to a CSV file
write.csv(data, "./data/GSE115469.Ensembl.csv", row.names = FALSE)

