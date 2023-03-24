rm(list = ls())

setwd("~/github/compbio_project/code")

cells_cohort_1_file = paste("../data/IntratumoralChangesAntiPD1/", 
                            "1863-counts_cells_cohort1.rds", sep="")
cohort1 <- readRDS(cells_cohort_1_file)
tcell_cohort1_file = paste("../data/IntratumoralChangesAntiPD1/", 
                            "1864-counts_tcell_cohort1.rds", sep="")
tcell_cohort1 <- readRDS(tcell_cohort1_file)
tcell_metadata_cohort1_file = paste("../data/IntratumoralChangesAntiPD1/", 
                           "1870-BIOKEY_metaData_tcells_cohort1_web.csv", sep="")
tcell_metadata_cohort1 <- read.csv(tcell_metadata_cohort1_file, sep="\t")
