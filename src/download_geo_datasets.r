if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)
gse <- getGEO("GSE81861")
gse_1 <- gse[[1]]

# exprs(gse[[1]])
# write.csv(exprs(gse[[1]]), file = "GSE81861_expression.csv")
# 
# metadata <- pData(gse[[1]])[1,]
