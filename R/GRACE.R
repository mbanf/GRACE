
#' internal function
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
build_estimate_rank_difference_cutoff <- function(mat.ranking=mat.ranking, tfs=tfs, select = "global", v.th.rank_cutoff = 0.05){
  
  v.rank_combs <- numeric()
  for(i in 1:length(tfs)){
    v.ranks <- mat.ranking[tfs[i],]
    v.ranks <- v.ranks[v.ranks > 0]
    if(length(v.ranks) > 1){  
      df.links <- as.data.frame(t(combn(v.ranks,2)), stringsAsFactors = FALSE)
      v.rank_combs <- c(v.rank_combs, abs(df.links$V1 - df.links$V2))
    }
  }
  v.cutoff <- as.numeric(quantile(v.rank_combs, v.th.rank_cutoff))
  return(v.cutoff) 
}


#' build link ranking based connectivity matrix
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param v.th.rank_cutoff default = 1
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
build_link_ranking_connectivity_matrix <- function(df.grn=df.grn, mat.RN=mat.RN, tfs=tfs, tgs=tgs, n.links=n.links, v.th.rank_cutoff=1, n.cores = 2){
  
  mat.RN.sparse <- Matrix(mat.RN)
  mat.ranking <- Matrix(0, nrow = length(tfs), ncol = n.links, dimnames = list(tfs, as.character(seq(1:n.links))))
  for(i in 1:length(tfs)){
    #cat("Processing... ", round(i/length(tfs) * 100, digits = 2) , "%", "\r"); flush.console()    
    df.grn.set <- subset(df.grn, df.grn$TF %in% tfs[i])
    rownames(df.grn.set) <- df.grn.set$link.id
    mat.ranking[i, df.grn.set$link.id] <- df.grn.set$v.grn  #df.grn.set$rank.global ##
  }
  
  v.links_tgs.map <- df.grn$link.id
  names(v.links_tgs.map) <- df.grn$tg.id
  v.links_tgs.map <- v.links_tgs.map[order(v.links_tgs.map)]
  
  v.cutoff <- build_estimate_rank_difference_cutoff(mat.ranking=mat.ranking, tfs=tfs, select = select, v.th.rank_cutoff = v.th.rank_cutoff)
  n.total <- length(tfs)
  n.each <- floor(n.total/n.cores)
  
  cl<-makeCluster(n.cores)
  registerDoParallel(cl)
  strt<-Sys.time()
  lst.adj <- foreach(r = 1:n.cores, .packages=c("Matrix", "reshape2")) %dopar% {    
    v.start <- n.each * (r - 1) + 1
    v.stop <- n.each * r
    tfs.set <- tfs[v.start:v.stop]
    mat.ranks <- Matrix(0, nrow = n.links, ncol = n.links, dimnames = list(as.character(seq(1:n.links)), as.character(seq(1:n.links))))
    for(i in 1:length(tfs.set)){
      cat("Processing... ", round(i/length(tfs.set) * 100, digits = 2) , "%", "\r"); flush.console()    
      v.ranks <- mat.ranking[tfs.set[i],]
      v.ranks <- v.ranks[v.ranks > 0]
      if(length(v.ranks) > 1){  
        df.links <- as.data.frame(t(combn(v.ranks,2)), stringsAsFactors = FALSE)
        v.rank_combs <- abs(df.links$V1 - df.links$V2) 
        idx <- which(v.rank_combs <= v.cutoff) #seq(1:length(v.rank_combs)) #
        if(length(idx) > 0){
          v.rank_combs <- v.rank_combs[idx]  
          df.links <- as.data.frame(t(combn(names(v.ranks),2)), stringsAsFactors = FALSE)
          df.links <- df.links[idx,]
          
          l1 <- as.numeric(df.links$V1)
          l2 <- as.numeric(df.links$V2)
          
          ind.pcc <- cbind(as.numeric(names(v.links_tgs.map)[l1]), as.numeric(names(v.links_tgs.map)[l2]))
          v.pccs <- mat.RN.sparse[ind.pcc]
          
          ind <- cbind(l1,l2)
          mat.ranks[ind] <- v.pccs / (1 + v.rank_combs)
          #mat.ranks[ind] <- 1 / (1 + v.rank_combs) # modified to test influence of Aranet
          #mat.ranks[ind] <- (v.rank_combs) / (1 +v.pccs) # modified to test inverse (negative test)
        }
      }
    }
    mat.ranks
  }
  stopCluster(cl)
  mat.ranks_links <- Reduce('+',lst.adj)
  
  ### residuals
  v.start <- n.each * n.cores + 1
  v.stop <- n.total 
  if(v.start <= v.stop){
    tfs.set <- tfs[v.start:v.stop]
    for(i in 1:length(tfs.set)){
      #cat("Processing... ", round(i/length(tfs.set) * 100, digits = 2) , "%", "\r"); flush.console()    
      v.ranks <- mat.ranking[tfs.set[i],]
      v.ranks <- v.ranks[v.ranks > 0]
      if(length(v.ranks) > 1){  
        df.links <- as.data.frame(t(combn(v.ranks,2)), stringsAsFactors = FALSE)
        v.rank_combs <- abs(df.links$V1 - df.links$V2) 
        idx <- which(v.rank_combs <= v.cutoff) # seq(1:length(v.rank_combs)) #
        if(length(idx) > 0){
          v.rank_combs <- v.rank_combs[idx]
          df.links <- as.data.frame(t(combn(names(v.ranks),2)), stringsAsFactors = FALSE)
          df.links <- df.links[idx,]
          l1 <- as.numeric(df.links$V1)
          l2 <- as.numeric(df.links$V2)
          
          ind.pcc <- cbind(as.numeric(names(v.links_tgs.map)[l1]), as.numeric(names(v.links_tgs.map)[l2]))
          v.pccs <- mat.RN.sparse[ind.pcc]
          
          ind <- cbind(l1,l2)
          mat.ranks_links[ind] <- v.pccs / (1 + v.rank_combs)
          #mat.ranks_links[ind] <- 1 / (1 + v.rank_combs) # modified to test influence of Aranet
          #mat.ranks_links[ind] <- (v.rank_combs) / (v.pccs + 1) # modified to test inverse (negative test)
        } 
      }
    }
  }
  mat.ranks_links.target_connectivity <- mat.ranks_links + t(mat.ranks_links)
  
  rm(mat.RN.sparse)
  rm(mat.ranking)
  rm(lst.adj)
  
  mat.ranks_links <- mat.ranks_links.target_connectivity
  
  return(mat.ranks_links)
}



#' construct meta network to propagate regulatory information
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
contruct_grace_meta_connectivity_network <- function(mat.grn=mat.grn,  mat.FunctionalAssociation=mat.FunctionalAssociation,
                                                     th.preselection = 0.0, n.cores = 2){
  
  mat.RN <- mat.FunctionalAssociation
  v.grn.significant <- th.preselection
  
  df.grn <- as.data.frame(as.table(mat.grn), stringsAsFactors = FALSE) 
  names(df.grn)[1:3] <- c("TF","TG","v.grn")
  df.grn <- subset(df.grn, df.grn$v.grn > 0)
  df.grn <- df.grn[order(-df.grn$v.grn),] 
  
  if(v.grn.significant < df.grn$v.grn[nrow(df.grn)]){
    v.grn.significant <- df.grn$v.grn[nrow(df.grn)]
  }
  df.grn <- subset(df.grn, df.grn$v.grn >= v.grn.significant)
  
  tfs <- unique(df.grn$TF) 
  tgs <- unique(df.grn$TG)
  
  v.grn.node <- df.grn$v.grn
  n.links <- nrow(df.grn)
  df.grn["link.id"] <- seq(1, n.links) 
  rownames(df.grn) <- df.grn$link.id
  
  tf.map <- new.env(hash=T, parent=emptyenv())
  for(j in 1:length(tfs)){ tf.map[[tfs[j]]] <- j}
  df.grn["tf.id"] <- 0
  for(j in 1:length(tfs)){
    #cat("Processing... ", round(j/length(tfs) * 100, digits = 2) , "%", "\r"); flush.console()   
    idx <- which(df.grn$TF == tfs[j]); df.grn$tf.id[idx] <- tf.map[[tfs[j]]]
  }
  
  tg.map <- new.env(hash=T, parent=emptyenv())
  for(j in 1:length(tgs)){ tg.map[[tgs[j]]] <- j}
  df.grn["tg.id"] <- 0
  for(j in 1:length(tgs)){    
    #cat("Processing... ", round(j/length(tgs) * 100, digits = 2) , "%", "\r"); flush.console()   
    idx <- which(df.grn$TG == tgs[j]); df.grn$tg.id[idx] <- tg.map[[tgs[j]]]
  }
  
  mat.link_id <- acast(df.grn, TF~TG, value.var="link.id")
  mat.link_id[is.na(mat.link_id)] <- 0; 
  mat.link_id <- as(mat.link_id, "CsparseMatrix")
  
  ### IMPORTANT - Reorder tfs and tgs in all relevant matrices
  mat.link_id <- mat.link_id[tfs, tgs]
  mat.RN <- mat.RN[tgs, tgs]; 
  
  v.links.grn <- df.grn$link.id
  
  adj.links <- build_link_ranking_connectivity_matrix(df.grn=df.grn, mat.RN=mat.RN, tfs=tfs, tgs=tgs, n.links=n.links, n.cores = n.cores, v.th.rank_cutoff=1)[v.links.grn, v.links.grn]
  rm(mat.RN)
  
  return(list(df.grn=df.grn, adj.links=adj.links, v.grn.significant=v.grn.significant))
}


#' extract meta modules
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
prepare_meta_modules <- function(df.grn=df.grn, adj.links=adj.links){
  
  mat.grn.specs <- as.matrix(df.grn)
  
  v.links.singletons <- rownames(adj.links)[which(rowSums(adj.links) == 0)] 

  ## remove singleton nodes
  adj.links.nonsingleton <- adj.links[!rownames(adj.links) %in% v.links.singletons, !colnames(adj.links) %in% v.links.singletons]
  #adj.links.singleton <- adj.links[rownames(adj.links) %in% v.links.singletons, colnames(adj.links) %in% v.links.singletons]
  
  ### STEP 2 - compute Markov clustering     
  #print("extract link modules...")
  g.coreg <- graph.adjacency(adj.links.nonsingleton, weighted = TRUE, mode = "undirected")
  clust.coreg <- clusters(g.coreg) # igraph based pre-MCL extraction of connected components s
  v.coreg_modules <- unique(clust.coreg$membership)
  n.crfs <- length(v.coreg_modules)
  
  lst.crfs <- vector(mode = "list", length = n.crfs)
  lst.crf.prim_based_links <- vector(mode = "list", length = n.crfs)
  lst.v.grn.node <- vector(mode = "list", length = n.crfs)
  lst.tf.ids <- vector(mode = "list", length = n.crfs)
  lst.tg.ids <- vector(mode = "list", length = n.crfs)
  
  for(i in 1:length(v.coreg_modules)){ 
    # cat("Processing... ", round(i/length(v.coreg_modules) * 100, digits = 2) , "%", "\r"); flush.console()  
    
    v.links <- names(which(clust.coreg$membership == v.coreg_modules[i]))
    names(v.links) <- seq(1:length(v.links))
    
    adj.links.set <- adj.links.nonsingleton[v.links,v.links]
    
    ru.adj.links <- adj.links.nonsingleton[v.links,v.links]
    ru.adj.links@x <- 1 - ru.adj.links@x
    g.ru <- graph.adjacency(ru.adj.links, weighted = TRUE, mode = "undirected")
    g.ru.mst <- minimum.spanning.tree(g.ru, algorithm= "prim")
    
    df.ru.edges <- get.edges(g.ru.mst, E(g.ru.mst))
    ru.adj.links <- Matrix(0, nrow = length(v.links), ncol = length(v.links), dimnames = list(v.links,v.links))
    
    l1 <- as.numeric(names(v.links[df.ru.edges[,1]]))
    l2 <- as.numeric(names(v.links[df.ru.edges[,2]]))
    ru.adj.links[cbind(l1,l2)] <- adj.links.set[cbind(l1,l2)]
    rm(adj.links.set)
    
    crf <- make.crf(ru.adj.links, n.states = 2)
    
    lst.crfs[[i]] <- crf
    lst.crf.prim_based_links[[i]] <- ru.adj.links[crf$edges]
    
    v.links <- names(which(clust.coreg$membership == v.coreg_modules[i]))
    mat.set <- mat.grn.specs[v.links,]
    v.grn.node <- as.numeric(mat.set[,3]) # instead of rank : 6
    lst.v.grn.node[[i]] <- v.grn.node
    
    tf.ids <- as.numeric(mat.set[,5])
    tg.ids <- as.numeric(mat.set[,6])
    
    lst.tf.ids[[i]] <- tf.ids
    lst.tg.ids[[i]] <- tg.ids
  } 
  
  return(list(v.coreg_modules=v.coreg_modules, lst.crfs=lst.crfs, lst.crf.prim_based_links=lst.crf.prim_based_links, g.coreg=g.coreg, 
              lst.v.grn.node=lst.v.grn.node, lst.tf.ids=lst.tf.ids, lst.tg.ids=lst.tg.ids))
}


#' generate bootstrap samples for analysis based on extracting sets of transcription factors
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
generate_bootstrapping_sample <- function(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = 0.632, b.validation = FALSE){
  
  mat.GS.bootstrap_sample <- truereg
  df.GS <- as.data.frame(which(mat.GS.bootstrap_sample == 1, arr.ind = TRUE))
  names(df.GS) <- c("TF", "TG")
  idx.tfs.gs <- unique(df.GS$TF)
  v.tfs.gs <- rownames(mat.GS.bootstrap_sample)[idx.tfs.gs]

  # sample from gold standard
  v.filter <- sample(nrow(df.GS), round(nrow(df.GS) * (n.sample_size)))
  
  if(b.validation){
    # select the complement for validation 
    v.filter.train <- v.filter
    v.filter <- seq(1:nrow(df.GS))[!seq(1:nrow(df.GS)) %in% v.filter.train]
  }
  
  # sample of the gold standard transcription factors (needs to be non-overlapping separated into training and testing)
  v.tfs.gs.bootstrap_sample <- unique(df.GS$TF[v.filter])
  v.tfs.gs.bootstrap_sample <- rownames(mat.GS.bootstrap_sample)[v.tfs.gs.bootstrap_sample]
  
  # sample from data (then replace ) - this can be identical to the training, as this is a semisupervised setting
  v.tfs <- (unique(lst.grace_support$df.grn$TF))
  names(v.tfs) <- unique(lst.grace_support$df.grn$tf.id)
  n.tfs <- length(v.tfs)
  v.tfs.ids <- unique(lst.grace_support$df.grn$tf.id)
  names(v.tfs.ids) <- (unique(lst.grace_support$df.grn$TF))
  
  n.tfs.bootstrap_sample <- round(n.tfs * n.sample_size, digits = 0)
  n.tfs.bootstrap_sample <- n.tfs.bootstrap_sample - length(v.tfs.gs.bootstrap_sample)

  v.tfs.bootstrap_sample.gs_plus <- sample(v.tfs[!v.tfs %in% v.tfs.gs], n.tfs.bootstrap_sample, replace = FALSE)
  n.tfs.bootstrap_sample.gs_plus <- length(v.tfs.bootstrap_sample.gs_plus)

  # >= 0.632 of features (regulators) - final selection
  v.tfs.bootstrap_sample <- unique(c(v.tfs.bootstrap_sample.gs_plus, v.tfs.gs.bootstrap_sample))
  v.tfs.filter <- v.tfs[!v.tfs %in% v.tfs.bootstrap_sample]
  
  mat.GS.bootstrap_sample[v.tfs.filter,] <- 0 # this might be re-designed
  
  # sample datasets using regulator bootstrap samples (adjust targets for max coreg bp pairs computation)
  lst.grace_support.bootstrap_sample <- lst.grace_support
  lst.grace_support.bootstrap_sample$df.grn <- subset(lst.grace_support.bootstrap_sample$df.grn, lst.grace_support.bootstrap_sample$df.grn$TF %in% v.tfs.bootstrap_sample)

  v.tgs.bootstrap_sample <- unique(lst.grace_support.bootstrap_sample$df.grn$TG)
  v.tgs.cofunctional_evidence <- colnames(mat.cofunctional_evidence)
  v.tgs.filter <- v.tgs.cofunctional_evidence[!v.tgs.cofunctional_evidence %in% v.tgs.bootstrap_sample]
  
  mat.cofunctional_evidence.bootstrap_sample <- mat.cofunctional_evidence
  mat.cofunctional_evidence.bootstrap_sample[v.tgs.filter,v.tgs.filter] <- 0

  ###

  v.link.bootstrap_sample <- lst.grace_support.bootstrap_sample$df.grn$link.id
  lst.grace_support.bootstrap_sample$adj.links <- lst.grace_support.bootstrap_sample$adj.links[v.link.bootstrap_sample,v.link.bootstrap_sample]
  
  ## more regulators than regulators?
  lst.modules.bootstrap_sample <- lst.modules
  idx.tfs_ids.bootstrap_sample <- unlist(lapply(lst.modules.bootstrap_sample$lst.tf.ids, function(m) unlist(m[[1]])))
  

  idx.tfs_ids.bootstrap_sample <- which(idx.tfs_ids.bootstrap_sample %in% as.numeric(v.tfs.ids[v.tfs.bootstrap_sample]))
  
  lst.modules.bootstrap_sample$v.coreg_modules <- lst.modules$v.coreg_modules[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.crfs <- lst.modules$lst.crfs[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.crf.prim_based_links <- lst.modules$lst.crf.prim_based_links[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.v.grn.node <- lst.modules$lst.v.grn.node[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.tf.ids <- lst.modules$lst.tf.ids[idx.tfs_ids.bootstrap_sample]
  lst.modules.bootstrap_sample$lst.tg.ids <- lst.modules$lst.tg.ids[idx.tfs_ids.bootstrap_sample]
  
  
  return(list(lst.grace_support.bootstrap_sample=lst.grace_support.bootstrap_sample, lst.modules.bootstrap_sample=lst.modules.bootstrap_sample, 
              lst.benchmark.bootstrap_sample=list(mat.GS.bootstrap_sample=mat.GS.bootstrap_sample, mat.cofunctional_evidence.bootstrap_sample=mat.cofunctional_evidence.bootstrap_sample)))
}





#' re-evaluate regulatory links of initial network using meta modules using Markov Random Fields
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
reevaluate_links_with_MarkovRandomFields <- function(lst.grace_support=lst.grace_support, lst.modules=lst.modules,
                                                     v.alpha=v.alpha, v.gamma=v.gamma, v.lambda=v.lambda,  b.prim = TRUE,
                                                     tfs = unique(lst.grace_support$df.grn$TF), tgs = unique(lst.grace_support$df.grn$TG)){
  
  th.gamma <- v.gamma
  th.alpha <- v.alpha
  
  df.grn <- lst.grace_support$df.grn
  adj.links <- lst.grace_support$adj.links
  th.preselection <- lst.grace_support$v.grn.significant 
  
  g.coreg <- lst.modules$g.coreg
  v.coreg_modules <- lst.modules$v.coreg_modules
  n.crfs <- length(v.coreg_modules)
  lst.crfs.prim <- lst.modules$lst.crfs
  lst.crf.prim_based_links <- lst.modules$lst.crf.prim_based_links
  lst.v.grn.node <- lst.modules$lst.v.grn.node
  
  lst.tf.ids <- lst.modules$lst.tf.ids
  lst.tg.ids <- lst.modules$lst.tg.ids
  
  #   mat.link_id <- acast(df.grn, TF~TG, value.var="link.id")
  #   mat.link_id[is.na(mat.link_id)] <- 0;  
  #   mat.link_id <- as(mat.link_id, "CsparseMatrix")
  #   
  #   ### IMPORTANT - Reorder tfs and tgs in all relevant matrices
  #   mat.link_id <- mat.link_id[tfs, tgs]
  #   
  # v.links.grn <- df.grn$link.id
  
  ## singleton nodes                         
  v.links.singletons <- rownames(adj.links)[which(rowSums(adj.links) == 0)] 
  
  ## remove singleton nodes
  adj.links.nonsingleton <- adj.links[!rownames(adj.links) %in% v.links.singletons, !colnames(adj.links) %in% v.links.singletons]

  clust.coreg <- clusters(g.coreg) # igraph based extraction of connected components s
  
  lst.crfs <- vector(mode = "list", length = n.crfs)
  for(i in 1:length(v.coreg_modules)){ 

    v.links <- names(which(clust.coreg$membership == v.coreg_modules[i]))
    adj.links.set <- adj.links.nonsingleton[v.links,v.links]
    
    v.links %in% rownames(adj.links.nonsingleton)
    v.links %in% rownames(adj.links)
    
    if(b.prim){
      crf <- lst.crfs.prim[[i]]
      v.edge_pots <- exp(v.lambda * lst.crf.prim_based_links[[i]]) # ru.adj.links[crf$edges])
    } else{
      crf <- make.crf(adj.links.set, n.states = 2)
      v.edge_pots <- exp(v.lambda * adj.links.set[crf$edges])
    }
    ## edge potentials
    for(j in 1:crf$n.edges){ # only lambda
      crf$edge.pot[[j]] <- matrix(c(v.edge_pots[j], 1, 1, v.edge_pots[j]), ncol = 2, byrow = FALSE)
    } 
    lst.crfs[[i]] <- crf
  }
  
  ### Compute MRF
  if(length(v.links.singletons) > 0){
    df.grn.singletons <- subset(df.grn, df.grn$link.id %in% v.links.singletons)
    #mat.p.crf_precomputed <- acast(df.grn.singletons, TF~TG, value.var = "rank.global")
    mat.p.crf_singleton <- acast(df.grn.singletons, TF~TG, value.var = "v.grn")
    mat.p.crf_singleton[is.na(mat.p.crf_singleton)] <- 0
    mat.p.crf_singleton <- exp(mat.p.crf_singleton - th.gamma)
    mat.p.crf_singleton[mat.p.crf_singleton <= 1] <- 0
    mat.p.crf_singleton[mat.p.crf_singleton > 1] <- 1
    #mat.p.crf_singleton <- as(mat.p.crf_singleton, "CsparseMatrix")
  }


  mat.p.crf <- matrix(0, nrow = length(tfs), ncol = length(tgs), dimnames = list(tfs, tgs)) 
  v.counter.crf <- numeric(length(v.coreg_modules))
  
  for(i in 1:length(v.coreg_modules)){ 
    v.grn.node <- lst.v.grn.node[[i]]
    crf <- lst.crfs[[i]]
    crf$node.pot[,1] <- exp(v.grn.node - th.gamma)
    crf$node.pot[,2] <- 1 # mean normalized to 0
    
    idx.gamma <- which(v.grn.node > th.alpha)
    #idx.gamma <- which(crf$node.pot[,1] > 1) # everything > 0.5 is automatically included (p = 1)
    clamped <- rep(0, length(crf$node.pot[,1]))
    clamped[idx.gamma] <- 1  
    
    if(!all(clamped == 1)){
      p.crf <- infer.conditional(crf, clamped, infer.tree)$node.bel[,1] # infer marginals / node beliefs 
    }else{
      v.counter.crf[i] <- 1
      p.crf <- rep(1.0, length(clamped))
    }
    
    p.crf <- ifelse(p.crf > 0.5, 1, 0) 
    tf.ids <- lst.tf.ids[[i]]
    tg.ids <- lst.tg.ids[[i]]
    ind = (tg.ids-1)*nrow(mat.p.crf) + tf.ids
    mat.p.crf[ind] <- p.crf
  
  }
  
  if(length(v.links.singletons) > 0){
    tf.pre <- rownames(mat.p.crf_singleton)
    tg.pre <- colnames(mat.p.crf_singleton)
    mat.p.crf[tf.pre, tg.pre] <- mat.p.crf[tf.pre, tg.pre] + mat.p.crf_singleton
  }
  
  if(all(v.counter.crf == 1)){
    print("Note: all links to be re-evaluated based on functional association data already in pre-selection, no changes")
  }

  return(list(mat.p.crf=mat.p.crf, th.gamma=th.gamma))
}



