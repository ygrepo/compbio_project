rm(list = ls())

# load
library(MuSiC)
library(SingleCellExperiment)

# Download Mouse bulk dataset from Github
Mouse.bulk.eset = readRDS('../data/Mousebulkeset.rds')

# Download Mouse single cell dataset from Github
Mousesub.sce = readRDS('../data/Mousesub_sce.rds')

levels(Mousesub.sce$cellType)

# Produce the first step information
Mousesub.basis = music_basis(Mousesub.sce, clusters = 'cellType', samples = 'sampleID', 
                             select.ct = c('Endo', 'Podo', 'PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC', 'Fib',
                                           'Macro', 'Neutro','B lymph', 'T lymph', 'NK'))

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(Mousesub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(Mousesub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')


clusters.type = list(C1 = 'Neutro', C2 = 'Podo', C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'), C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))

cl.type = as.character(Mousesub.sce$cellType)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
Mousesub.sce$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
s.mouse = unlist(clusters.type)
s.mouse


load('../data/IEmarkers.RData')
# This RData file provides two vectors of gene names Epith.marker and Immune.marker

# We now construct the list of group marker
IEmarkers = list(NULL, NULL, Epith.marker, Immune.marker)
names(IEmarkers) = c('C1', 'C2', 'C3', 'C4')
# The name of group markers should be the same as the cluster names

Est.mouse.bulk = music_prop.cluster(bulk.mtx = exprs(Mouse.bulk.eset), sc.sce = Mousesub.sce, group.markers = IEmarkers, clusters = 'cellType', groups = 'clusterType', samples = 'sampleID', clusters.type = clusters.type)

