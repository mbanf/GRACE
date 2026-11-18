
#' load example datasets for A. thaliana development network inference and evaluation
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
request_dataset <- function(species = "arabidopsis", method = "GENIE3"){
  
  lst.data <- vector(mode = "list", length = 8)
  if(species == "arabidopsis"){
    
    print("loading A. thaliana datasets...")
    
    print("ATRM & AGRIS binding evidence...")
    df.AtRegNet.novel <- read.table("athaliana_datasets/reg_net.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    df.AtRegNet.novel <- subset(df.AtRegNet.novel, df.AtRegNet.novel$V8 == "Confirmed") # AtRegNet gold Standard
    df.AtRegNet.novel <- df.AtRegNet.novel[,c(2,5)]
    names(df.AtRegNet.novel) <- c("TF", "Target")
    df.AtRegNet.novel$TF <- toupper(df.AtRegNet.novel$TF)
    df.AtRegNet.novel$Target <- toupper(df.AtRegNet.novel$Target)
    
    mat.AtRegNet <- acast(df.AtRegNet.novel, TF~Target)
    mat.AtRegNet[is.na(mat.AtRegNet)] <- 0
    mat.AtRegNet[mat.AtRegNet != 0] <- 1
    class(mat.AtRegNet) <- "numeric"   
    
    mat.AtRegNet <- mat.AtRegNet[names(which(rowSums(mat.AtRegNet) < 2000)),] # remove 
    
    
    df.ATRM <- read.table("athaliana_datasets/ATRM.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    mat.ATRM <- acast(df.ATRM, TF~Target)
    mat.ATRM[is.na(mat.ATRM)] <- 0
    mat.ATRM[mat.ATRM != 0] <- 1
    class(mat.ATRM) <- "numeric"   
    mat.GS <- mat.ATRM
    
    # the filtered regulatory network based on integrating binding with developmental gene expression ( > 0.95%le )
    mat.conserved_motifs_GENIE3 <- readRDS("athaliana_datasets/mat.conserved_motifs_plus_genie3_development.rds")
    v.tfs <- rownames(mat.conserved_motifs_GENIE3)
    v.tgs <- colnames(mat.conserved_motifs_GENIE3)
    v.gns <- unique(c(v.tfs, v.tgs))
    
    print("Aranet...")
    # benchmark used to evaluate Aranet V2 (http://www.inetbio.org/aranet/), experimental gene ontology annotations
    df.aranet_benchmark <- read.table("athaliana_datasets/AraNet.v2.goldstandard.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(df.aranet_benchmark) <- c("gn.a", "gn.b")
    df.aranet_benchmark <- subset(df.aranet_benchmark, df.aranet_benchmark$gn.a %in% v.gns)
    df.aranet_benchmark <- subset(df.aranet_benchmark, df.aranet_benchmark$gn.b %in% v.gns)
    df.aranet_benchmark["val"] <- 1
    mat.aranet_benchmark <- acast(df.aranet_benchmark, gn.a~gn.b, value.var = "val")
    mat.aranet_benchmark[is.na(mat.aranet_benchmark)] <- 0
    class(mat.aranet_benchmark) <- "numeric"   
    mat.aranet_benchmark <- as(mat.aranet_benchmark, "CsparseMatrix")
    
    mat.grn.aranet_benchmark <- Matrix(0, nrow = length(v.tfs), ncol = length(v.tgs), dimnames = list(v.tfs, v.tgs))
    gns.rows <- intersect(rownames(mat.grn.aranet_benchmark), rownames(mat.aranet_benchmark))
    gns.cols <- intersect(colnames(mat.grn.aranet_benchmark), colnames(mat.aranet_benchmark))
    mat.grn.aranet_benchmark[gns.rows, gns.cols] <- mat.aranet_benchmark[gns.rows, gns.cols]
    
    mat.tgs.aranet_benchmark <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
    gns.rows <- intersect(rownames(mat.tgs.aranet_benchmark), rownames(mat.aranet_benchmark))
    gns.cols <- intersect(colnames(mat.tgs.aranet_benchmark), colnames(mat.aranet_benchmark))
    mat.tgs.aranet_benchmark[gns.rows, gns.cols] <- mat.aranet_benchmark[gns.rows, gns.cols]
    
    mat.tfs.aranet_benchmark <- Matrix(0, nrow = length(v.tfs), ncol = length(v.tfs), dimnames = list(v.tfs, v.tfs))
    gns.rows <- intersect(rownames(mat.tfs.aranet_benchmark), rownames(mat.aranet_benchmark))
    gns.cols <- intersect(colnames(mat.tfs.aranet_benchmark), colnames(mat.aranet_benchmark))
    mat.tfs.aranet_benchmark[gns.rows, gns.cols] <- mat.aranet_benchmark[gns.rows, gns.cols]
    
    
    # most recent release of Aranet V2 (http://www.inetbio.org/aranet/)
    df.aranet <- read.table("athaliana_datasets/AraNet.v2.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(df.aranet) <- c("gn.a", "gn.b", "value")
    df.aranet <- subset(df.aranet, df.aranet$gn.a %in% v.tgs)
    df.aranet <- subset(df.aranet, df.aranet$gn.b %in% v.tgs)
    mat.aranet <- acast(df.aranet, gn.a~gn.b, value.var = "value")
    mat.aranet[is.na(mat.aranet)] <- 0
    class(mat.aranet) <- "numeric"  
    mat.aranet <- as(mat.aranet, "CsparseMatrix")
    
    # expand aranet matrix to cover all tgs
    mat.RN <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
    gns.rows <- intersect(rownames(mat.RN), rownames(mat.aranet))
    gns.cols <- intersect(colnames(mat.RN), colnames(mat.aranet))
    mat.RN[gns.rows, gns.cols] <- mat.aranet[gns.rows, gns.cols]
    
    # update here to account for symmetry
    
    mat.RN <- mapping_intervals(mat.RN, min(mat.RN), max(mat.RN), 0, 1) 
    
    
    ### load additional validation data
    
    
    # Subcellular co-localization
    df.suba <- read.csv("athaliana_datasets/Suba3-2015-2-10_7-18.csv", header = TRUE, stringsAsFactors = FALSE)
    df.suba$AGI.identifier <- gsub("\\..*", "", df.suba$AGI.identifier)
    df.suba <- subset(df.suba, df.suba$AGI.identifier %in% v.gns)
    df.suba <- df.suba[,c("AGI.identifier", "location_gfp", "location_swissprot", "location_ms.ms", "location_tair")] # subset to experimental data only 
    vec.compartments.unique <- character()
    
    for(k in 1:nrow(df.suba)){
      cat("Processing... ", round(k/nrow(df.suba) * 100, digits = 2) , "%", "\r"); flush.console()  
      gn.id <- df.suba$AGI.identifier[k]
      if(gn.id %in% v.gns){
        
        vec.compartments.unique <- c(vec.compartments.unique, as.character(unique(c(unlist(strsplit(df.suba$location_gfp[k], ";")),
                                                                                    unlist(strsplit(df.suba$location_swissprot[k], ";")),
                                                                                    unlist(strsplit(df.suba$location_ms.ms[k], ";")),
                                                                                    unlist(strsplit(df.suba$location_tair[k], ";"))))))
      }
    }
    
    vec.compartments.unique <- sapply(vec.compartments.unique,  function(m) gsub("\\:.*","", m))
    vec.compartments.unique <- unique(vec.compartments.unique)
    
    mat.subacon <- Matrix(0, nrow = length(v.tgs), ncol = length(vec.compartments.unique), dimnames = list(v.tgs, vec.compartments.unique))
    for(k in 1:nrow(df.suba)){
      cat("Processing... ", round(k/nrow(df.suba) * 100, digits = 2) , "%", "\r"); flush.console()  
      gn.id <- df.suba$AGI.identifier[k]
      if(gn.id %in% v.tgs){
        v.compartments.k <- as.character(unique(c(unlist(strsplit(df.suba$location_gfp[k], ";")),
                                                  unlist(strsplit(df.suba$location_swissprot[k], ";")),
                                                  unlist(strsplit(df.suba$location_ms.ms[k], ";")),
                                                  unlist(strsplit(df.suba$location_tair[k], ";")))))
        v.compartments.k <- as.character(sapply(v.compartments.k,  function(m) gsub("\\:.*","", m)))
        if(length(v.compartments.k) > 0){
          mat.subacon[gn.id, v.compartments.k] <- 1    
        }
      }
    }
    
    # generate co-localization matrix (between target genes)
    mat.subacon <- mat.subacon[,c("mitochondrion", "nucleus", "plastid", "peroxisome","extracellular","endoplasmic reticulum", "vacuole", "golgi")]
    mat.subacon <- as.matrix(mat.subacon)
    mat.MR <- jaccard_count(as.matrix(mat.subacon))
    mat.MR <- as.matrix(mat.MR); diag(mat.MR) <- 0;
    rownames(mat.MR) <- colnames(mat.MR) <- rownames(mat.subacon) 
    mat.MR <- as(mat.MR, "CsparseMatrix")
    
    mat.MR@x[mat.MR@x < 1] <- 0
    mat.MR@x[mat.MR@x >= 1] <- 1
    
    mat.MR.subacon.tgs <- mat.MR
    
    
    ## ARACYC - same metabolic pathway cooccurence
    # Optional: validation-only dataset (not used for model training or the ensemble).
    # If the file is absent, an all-zero matrix is used and the AraCyc validation column reads 0.
    file.aracyc <- "athaliana_datasets/aracyc-12.0.1-pathways.20150407"
    mat.aracyc.tgs <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
    v.pwys <- character(0)
    if(file.exists(file.aracyc)){
      df.aracyc <- read.table(file.aracyc, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
      df.aracyc <- subset(df.aracyc, df.aracyc$Gene.id != "unknown")
      v.pwys <- unique(df.aracyc$Pathway.id)
    }else{
      warning("AraCyc pathway file not found - skipping metabolic pathway validation: ", file.aracyc)
    }
    for(p in seq_along(v.pwys)){
      v.gns.p <- subset(df.aracyc, df.aracyc$Pathway.id == v.pwys[p])$Gene.id    
      v.gns.p <- intersect(v.gns.p, v.tgs)
      if(length(v.gns.p) >= 2){
        mat.aracyc.tgs[v.gns.p,v.gns.p] <- 1  
      }
    }
    
    return(list(mat.conserved_motif_plus_expression=mat.conserved_motifs_GENIE3[v.tfs,], 
                v.tfs=v.tfs, v.tgs=v.tgs, mat.GS = mat.GS, mat.FunctionalAssociation=mat.RN,
                mat.GO.tgs=mat.tgs.aranet_benchmark, 
                mat.MR.subacon.tgs=mat.MR.subacon.tgs, 
                mat.AtRegNet=mat.AtRegNet, 
                mat.aracyc.tgs=mat.aracyc.tgs))
    
  }else if(species == "drosophila"){
    
    print("loading D. melanogaster datasets...")
    
    print("Redfly binding evidence...")
  
    df.GS <- read.table("drosophila_datasets/redfly.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    names(df.GS) <- c("TF", "Target", "V3")
    mat.GS <- acast(df.GS, TF~Target, value.var = "V3")
    mat.GS[is.na(mat.GS)] <- 0
    mat.GS[mat.GS != 0] <- 1
    class(mat.GS) <- "numeric"   
    
    
    #     df.GS <- read.table("drosophila_datasets/drosophila_network.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    #     names(df.GS) <- c("TF", "Target", "V3")
    #     mat.GS <- acast(df.GS, TF~Target, value.var = "V3")
    #     mat.GS[is.na(mat.GS)] <- 0
    #     mat.GS[mat.GS != 0] <- 1
    #     class(mat.GS) <- "numeric"   
    #     
    if(method == "GENIE3"){
      mat.conserved_motif_plus_expression <- readRDS("drosophila_datasets/mat.grn.rapid.drosophila.rds")  
    }else if(method == "CLR"){
      mat.conserved_motif_plus_expression <- readRDS("drosophila_datasets/mat.grn.clr_drosophila.rds")  
    }else if(method == "ARACNE"){
      mat.conserved_motif_plus_expression <- readRDS("drosophila_datasets/mat.grn.aracne_drosophila.rds")  
    }
    
    v.tfs <- rownames(mat.conserved_motif_plus_expression)
    v.tgs <- colnames(mat.conserved_motif_plus_expression)
    v.gns <- unique(c(v.tfs, v.tgs))
  
    # the filtered regulatory network based on integrating binding with developmental gene expression ( > 0.95%le )
    
    #     df.grn <- read.table("drosophila_datasets/flynet_unsupervised_0.02.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    #     names(df.grn) <- c("TF", "Target", "v.grn")
    #     mat.grn <- acast(df.grn, TF~Target, value.var = "v.grn")
    #     mat.grn[is.na(mat.grn)] <- 0
    #     #mat.grn[mat.grn != 0] <- 1
    #     class(mat.grn) <- "numeric"
    #   
    #     v.tfs <- rownames(mat.grn)
    #     v.tgs <- colnames(mat.grn)
    #     v.gns <- unique(c(v.tfs, v.tgs))
    
    print("Flynet...")
    df.flynet_benchmark <- read.table("drosophila_datasets/FlyNet.GS.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(df.flynet_benchmark) <- c("gn.a", "gn.b")
    df.flynet_benchmark <- subset(df.flynet_benchmark, df.flynet_benchmark$gn.a %in% v.gns)
    df.flynet_benchmark <- subset(df.flynet_benchmark, df.flynet_benchmark$gn.b %in% v.gns)
    df.flynet_benchmark["val"] <- 1
    mat.flynet_benchmark <- acast(df.flynet_benchmark, gn.a~gn.b, value.var = "val")
    mat.flynet_benchmark[is.na(mat.flynet_benchmark)] <- 0
    class(mat.flynet_benchmark) <- "numeric"   
    mat.flynet_benchmark <- as(mat.flynet_benchmark, "CsparseMatrix")
 
    mat.tgs.flynet_benchmark <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
    gns.rows <- intersect(rownames(mat.tgs.flynet_benchmark), rownames(mat.flynet_benchmark))
    gns.cols <- intersect(colnames(mat.tgs.flynet_benchmark), colnames(mat.flynet_benchmark))
    mat.tgs.flynet_benchmark[gns.rows, gns.cols] <- mat.flynet_benchmark[gns.rows, gns.cols]
    
    #### 
    
    # most recent release of Aranet V2 (http://www.inetbio.org/aranet/)
    df.flynet <- read.table("drosophila_datasets/FlyNet.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(df.flynet) <- c("gn.a", "gn.b", "value")
    df.flynet <- subset(df.flynet, df.flynet$gn.a %in% v.tgs)
    df.flynet <- subset(df.flynet, df.flynet$gn.b %in% v.tgs)
    mat.flynet <- acast(df.flynet, gn.a~gn.b, value.var = "value")
    mat.flynet[is.na(mat.flynet)] <- 0
    class(mat.flynet) <- "numeric"  
    mat.flynet <- as(mat.flynet, "CsparseMatrix")
    
    # expand aranet matrix to cover all tgs
    mat.RN <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
    gns.rows <- intersect(rownames(mat.RN), rownames(mat.flynet))
    gns.cols <- intersect(colnames(mat.RN), colnames(mat.flynet))
    mat.RN[gns.rows, gns.cols] <- mat.flynet[gns.rows, gns.cols]
    
    mat.RN <- mapping_intervals(mat.RN, min(mat.RN), max(mat.RN), 0, 1) 
    
    
    # Hi-C
    df.HiC <- read.table("drosophila_datasets/long_range_interactions.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(df.HiC) <- c("gn.a", "gn.b", "value")
    df.HiC <- subset(df.HiC, df.HiC$gn.a %in% v.tgs)
    df.HiC <- subset(df.HiC, df.HiC$gn.b %in% v.tgs)
    mat.HiC <- acast(df.HiC, gn.a~gn.b, value.var = "value")
    mat.HiC[is.na(mat.HiC)] <- 0
    mat.HiC[mat.HiC != 0] <- 1
    class(mat.HiC) <- "numeric"  
    mat.HiC.complete <- as(mat.HiC, "CsparseMatrix")
    
    mat.HiC <- Matrix(0, nrow = length(v.tgs), ncol = length(v.tgs), dimnames = list(v.tgs, v.tgs))
    gns.rows <- intersect(rownames(mat.HiC), rownames(mat.HiC.complete))
    gns.cols <- intersect(colnames(mat.HiC), colnames(mat.HiC.complete))
    mat.HiC[gns.rows, gns.cols] <- mat.HiC.complete[gns.rows, gns.cols]
    
    mat.HiC <- mapping_intervals(mat.HiC, min(mat.HiC), max(mat.HiC), 0, 1) 
    
    
    # use Chip evidence 
    df.Chip <- read.table("drosophila_datasets/chip_net.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
    names(df.Chip) <- c("gn.a", "gn.b", "value")
    df.Chip <- subset(df.Chip, df.Chip$gn.a %in% v.tfs)
    df.Chip <- subset(df.Chip, df.Chip$gn.b %in% v.tgs)
    mat.Chip <- acast(df.Chip, gn.a~gn.b, value.var = "value")
    mat.Chip[is.na(mat.Chip)] <- 0
    mat.Chip[mat.Chip != 0] <- 1
    class(mat.Chip) <- "numeric"  
    mat.Chip.complete <- as(mat.Chip, "CsparseMatrix")
    
    mat.Chip <- Matrix(0, nrow = length(v.tfs), ncol = length(v.tgs), dimnames = list(v.tfs, v.tgs))
    gns.rows <- intersect(rownames(mat.Chip), rownames(mat.Chip.complete))
    gns.cols <- intersect(colnames(mat.Chip), colnames(mat.Chip.complete))
    mat.Chip[gns.rows, gns.cols] <- mat.Chip.complete[gns.rows, gns.cols]
    
    mat.Chip <- mapping_intervals(mat.Chip, min(mat.Chip), max(mat.Chip), 0, 1) 
    

    
    # tissue location - GO 
    df.imago <- read.table("drosophila_datasets/imago.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    
    vec.compartments.unique <- unique(df.imago$V2)
    
    mat.tissue_localization <- Matrix(0, nrow = length(v.tgs), ncol = length(vec.compartments.unique), dimnames = list(v.tgs, vec.compartments.unique))
    for(k in 1:nrow(df.imago)){
      cat("Processing... ", round(k/nrow(df.imago) * 100, digits = 2) , "%", "\r"); flush.console()  
      gn.id <- df.imago$V1[k]
      if(gn.id %in% v.tgs){
        mat.tissue_localization[gn.id, df.imago$V2[k]] <- 1    
      }
    }
    
    # generate co-localization matrix (between target genes)
    mat.tissue_localization <- as.matrix(mat.tissue_localization)
    mat.MR <- jaccard_count(as.matrix(mat.tissue_localization))
    mat.MR <- as.matrix(mat.MR); diag(mat.MR) <- 0;
    rownames(mat.MR) <- colnames(mat.MR) <- rownames(mat.tissue_localization) 
    mat.MR <- as(mat.MR, "CsparseMatrix")
    
    mat.MR@x[mat.MR@x < 1] <- 0
    mat.MR@x[mat.MR@x >= 1] <- 1
    
    mat.tissue_localization.tgs <- mat.MR
    
    
    # redfly network 
    #     df.GS <- read.table("drosophila_datasets/redfly.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    #     names(df.GS) <- c("TF", "Target", "V3")
    #     mat.GS <- acast(df.GS, TF~Target, value.var = "V3")
    #     mat.GS[is.na(mat.GS)] <- 0
    #     mat.GS[mat.GS != 0] <- 1
    #     class(mat.GS) <- "numeric"   
    #     
    mat.REDfly<- Matrix(0, nrow = length(v.tfs), ncol = length(v.tgs), dimnames = list(v.tfs, v.tgs))
    gns.rows <- intersect(rownames(mat.REDfly), rownames(mat.GS))
    gns.cols <- intersect(colnames(mat.REDfly), colnames(mat.GS))
    mat.REDfly[gns.rows, gns.cols] <- mat.GS[gns.rows, gns.cols]
    #     
    

    return(list(mat.conserved_motif_plus_expression=mat.conserved_motif_plus_expression[v.tfs,], 
                v.tfs=v.tfs, v.tgs=v.tgs, mat.GS = mat.REDfly, mat.FunctionalAssociation=mat.RN,
                mat.GO.tgs=mat.tgs.flynet_benchmark,
                mat.HiC=mat.HiC,
                mat.Chip=mat.Chip,
                mat.tissue_localization.tgs=mat.tissue_localization.tgs
                #mat.REDfly=mat.REDfly
                ))  
    
  }
  

}
