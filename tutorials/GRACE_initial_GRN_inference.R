
# Contact: Michael Banf, michael.banf@gmx.net	
# 
# This is the documentation for the GRACE algorithm. The implementation is a research prototype and is provided “as is”. No warranties or guarantees of any kind are provided. Do not distribute the GRACE R code or use it other than for your own research without permission by the author. 
# 
# GRACE has been written in R Version 3.2.2 - you might have to update to this version here https://cran.r-project.org



rm(list=ls()) # clear workspace 

# set path to the GRACE directory
setwd("/home/mbanf/Documents/Computational_Biology/Projects/GRACE_release/")

# general functions
source("GRACE_helperfunctions.R")

cat("\014")  
message("request and load libraries..")
request_libraries()

th.grnCutoff <- 0.999 # initial cutoff for random forest regression on gene expression network
n.cpus <- 2 # number of available cpus for parallel random forest regression

m.expression <- ... # load expression matrix for grn inference - matrix - rows (genes) x cols (conditions)

# load set of transcription factors (here Arabidopsis) to be selected as regulators for inference
v.tfs <- ... # character vector of transcription factor genes

# compute random forest regression based grn 
m.rf <- compute_randomforest_based_GRN(m.expression, K="sqrt", nb.trees=1000, tfs=v.tfs, tgs = rownames(m.expression), importance.measure="IncNodePurity", seed=1234, n.cores = n.cpus)
# store network from gene expression 
saveRDS(m.rf, "tmp/m.grn_randomforest.rds")

m.rf <- readRDS("tmp/m.grn_randomforest.rds")
th.rf <- quantile(m.rf, th.grnCutoff) # keep only top predictions according to th.grnCutoff
m.rf[m.rf < th.rf] <- 0


df.dna_binding <- read.table(..., header  = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1:2] # load DNA binding set (dataframe with 2 columns)
names(df.dna_binding) <- c("gn.a", "gn.b")
m.dna_binding <- acast(df.dna_binding, gn.a~gn.b)
m.dna_binding[is.na(m.dna_binding)] <- 0
m.dna_binding[m.dna_binding != 0] <- 1
class(m.dna_binding) <- "numeric"  


v.tgs <- intersect(rownames(m.expression), colnames(m.dna_binding))
v.tfs <- intersect(v.tfs, rownames(m.dna_binding))


m.initial_grn <- m.rf[v.tfs, v.tgs] * m.dna_binding[v.tfs, v.tgs] # filter top predictions based on DNA binding
saveRDS(m.initial_grn, "m.initial_grn.rds")