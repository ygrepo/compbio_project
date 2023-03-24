rm(list = ls())
#BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)

listEnsembl()
ensembl <- useMart('ensembl')
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)



gpl.file <- 'GSE*.Ensembl.csv'
gpl.table <- read.csv(gpl.file,row.names=1,header=T)

gene.list <- row.names(gpl.table)


# test
gene_to_genome <- getBM(values = gene.list,
                 filters = 'ensembl_gene_id',
                 mart = ensembl,
                 attributes = c('ensembl_gene_id','external_gene_name'))

write.csv(gene_to_genome,'GSE*.Annotation.csv',row.names=FALSE)


