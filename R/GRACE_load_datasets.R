request_dataset <- function(){

  # set paths to novel datasets 
  lst.data <- vector(mode = "list", length = 8)

  message("loading datasets")
  
  # load regulatory evidence (2 column dataframe with experimental regulatory links)
  df.GS <- read.table(..., header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  m.GS <- acast(df.GS, TF~Target)
  m.GS[is.na(m.GS)] <- 0
  m.GS[m.GS != 0] <- 1
  class(m.GS) <- "numeric"   
  
  ### load DNA binding set 
  
  # the filtered regulatory network based on integrating binding with developmental gene expression ( as computed in GRACE_initieal_GRN_inference.R )
  m.initial_grn <- readRDS("m.initial_grn.rds")
  v.tfs <- rownames(m.initial_grn)
  v.tgs <- colnames(m.initial_grn)
  v.gns <- unique(c(v.tfs, v.tgs))
  
  # load function evidence  (2 column dataframe with experimental gene ontology annotations)
  df.cofunctional_benchmark <- read.table(..., header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df.cofunctional_benchmark) <- c("gn.a", "gn.b")
  df.cofunctional_benchmark <- subset(df.cofunctional_benchmark, df.cofunctional_benchmark$gn.a %in% v.gns)
  df.cofunctional_benchmark <- subset(df.cofunctional_benchmark, df.cofunctional_benchmark$gn.b %in% v.gns)
  df.cofunctional_benchmark["val"] <- 1
  m.cofunctional_benchmark <- acast(df.cofunctional_benchmark, gn.a~gn.b, value.var = "val")
  m.cofunctional_benchmark[is.na(m.cofunctional_benchmark)] <- 0
  class(m.cofunctional_benchmark) <- "numeric"   
  m.cofunctional_benchmark <- as(m.cofunctional_benchmark, "CsparseMatrix")
  
  mat.tgs.cofunctional_benchmark <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
  gns.rows <- intersect(rownames(mat.tgs.cofunctional_benchmark), rownames(m.cofunctional_benchmark))
  gns.cols <- intersect(colnames(mat.tgs.cofunctional_benchmark), colnames(m.cofunctional_benchmark))
  mat.tgs.cofunctional_benchmark[gns.rows, gns.cols] <- m.cofunctional_benchmark[gns.rows, gns.cols]
  
  # load cofunctional network (2 column dataframe)
  df.cofunctional <- read.table(..., header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df.cofunctional) <- c("gn.a", "gn.b", "value")
  df.cofunctional <- subset(df.cofunctional, df.cofunctional$gn.a %in% v.tgs)
  df.cofunctional <- subset(df.cofunctional, df.cofunctional$gn.b %in% v.tgs)
  m.cofunctional <- acast(df.cofunctional, gn.a~gn.b, value.var = "value")
  m.cofunctional[is.na(m.cofunctional)] <- 0
  class(m.cofunctional) <- "numeric"  
  m.cofunctional <- as(m.cofunctional, "CsparseMatrix")
  
  # expand aranet matrix to cover all tgs
  mat.RN <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
  gns.rows <- intersect(rownames(mat.RN), rownames(m.cofunctional))
  gns.cols <- intersect(colnames(mat.RN), colnames(m.cofunctional))
  mat.RN[gns.rows, gns.cols] <- m.cofunctional[gns.rows, gns.cols]
  
  # update here to account for symmetry
  
  mat.RN <- mapping_intervals(mat.RN, min(mat.RN), max(mat.RN), 0, 1) 
  
  
  return(list(mat.conserved_motif_plus_expression=m.initial_grn[v.tfs,], 
              v.tfs=v.tfs, v.tgs=v.tgs, mat.GS = m.GS, mat.FunctionalAssociation=mat.RN,
              mat.GO.tgs=mat.tgs.cofunctional_benchmark))
  
}
