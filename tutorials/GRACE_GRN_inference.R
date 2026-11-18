
#setwd("/Users/michaelbanf/Documents/postdoctoral_work/Projects/GRACE_release/")
setwd("/home/mbanf/Documents/Computational_Biology/Projects/GRACE_release/")


library(reshape2)


#install.packages("iRafNet")
# library(iRafNet)
library(glasso)
library(GeneNet)
library(predictionet)

# most recent release of Aranet V2 (http://www.inetbio.org/aranet/)
if(FALSE){
  
  df.aranet <- read.table("athaliana_datasets/AraNet.v2.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df.aranet) <- c("gn.a", "gn.b", "value")
  mat.aranet <- acast(df.aranet, gn.a~gn.b, value.var = "value")
  mat.aranet[is.na(mat.aranet)] <- 0
  class(mat.aranet) <- "numeric"  
  #mat.aranet <- as(mat.aranet, "CsparseMatrix")
  
  m.expression <- readRDS("mat.expression_development_arabidopsis.rds")
  m.expression <- t(m.expression)
  
  mat.cns_binding <- readRDS("mat.cns_binding.rds")
  
  v.tfs <- read.table("athaliana_datasets/Ath_tf_ptfb.txt", header  = TRUE, sep = "\t", stringsAsFactors = FALSE)
  v.tfs <- unique(gsub("\\..*", "", v.tfs$gene_model))
  
}else{

  df.aranet <- read.table("drosophila_datasets/FlyNet.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df.aranet) <- c("gn.a", "gn.b", "value")
  mat.aranet <- acast(df.aranet, gn.a~gn.b, value.var = "value")
  mat.aranet[is.na(mat.aranet)] <- 0
  class(mat.aranet) <- "numeric"  
  #mat.aranet <- as(mat.aranet, "CsparseMatrix")

  # the filtered regulatory network based on integrating binding with developmental gene expression ( > 0.95%le )
  df.gene_expression <- read.table("drosophila_datasets/microarray.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  m.expression <- data.matrix(df.gene_expression[,2:ncol(df.gene_expression)], rownames.force = NA)
  rownames(m.expression) <- df.gene_expression[,1]
  m.expression <- t(m.expression)
  
  
  library(parmigene)
  m.mi <- parmigene::knnmi.all(t(m.expression))
  m.clr <- parmigene::clr(m.mi)
  saveRDS(m.clr, "drosophila_datasets/clr_1120.rds")
  
  #   library(parmigene)
  #   m.mi <- parmigene::knnmi.all(m.expression)
  #   m.clr <- parmigene::clr(m.mi)
  #   saveRDS(m.clr, "athaliana_datasets/clr.rds")
  #   saveRDS(m.clr[v.tfs,], "athaliana_datasets/mat.clr_arabidopsis.rds")
  #   
  df.conserved_motifs <- read.table("motif_regnet_highestconf.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df.conserved_motifs) <- c("gn.a", "gn.b", "value")
  mat.cns_binding <- acast(df.conserved_motifs, gn.a~gn.b, value.var = "value")
  mat.cns_binding[is.na(mat.cns_binding)] <- 0
  mat.cns_binding[mat.cns_binding != 0] <- 1
  class(mat.cns_binding) <- "numeric"  
  
  
  v.tfs <- read.table("drosophila_datasets/tf.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
  
  
}
  


