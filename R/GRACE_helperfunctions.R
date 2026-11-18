
plot_fmeasure_scores <- function(df.rate_density){
  
  max.gs.grn <- max(df.rate_density$n.gs.grn)
  max.coreg.tgs.bp <- max(df.rate_density$n.coreg.bp_tg_pairs)
  max(df.rate_density$n.tg_pairs.bp)
  
  # extract ratio of precision x recall
  Precision.gs.grn <- (df.rate_density$n.gs.grn/df.rate_density$n.pred)
  #Recall.coreg.tgs.bp <- df.rate_density$n.coreg.tg_pairs.BP/max(df.rate_density$n.coreg.tg_pairs.BP)
  
  ratio.gs.grn_vs_bp.normalization <- 1 / max(Precision.gs.grn) #/ max(Recall.coreg.tgs.bp)
  max.coreg.bp <- max(df.rate_density$n.coreg.bp_tg_pairs)
  
  ####
  x <- df.rate_density$n.pred / max(df.rate_density$n.pred)
  
  y1 <- (df.rate_density$n.gs.grn/df.rate_density$n.pred)
  y1[is.na(y1)] <- 0
  y1 <- y1/max(y1) # normalization
  #plot(x,y1, type = "l") # training (ATRM)
  
  y2 <- df.rate_density$n.coreg.bp_tg_pairs/max(df.rate_density$n.coreg.bp_tg_pairs)
  y2[is.na(y2)] <- 0
  y2 <- y2#/max(y2)
  #lines(x, y2, col = "red") # training (GO-BP TG<->TG)
  
  beta.precision = 0.5
  fscore_beta <- (1 + beta.precision^2)*(y1*y2/((beta.precision^2 * y1)+y2))
  fscore_beta[is.na(fscore_beta)] <- 0
  fscore_beta_05 <- fscore_beta#/max(fscore_beta)
  
  beta.traditional = 1 #(44/10087)/(175/98596)
  fscore_beta <- (1 + beta.traditional^2)*(y1*y2/((beta.traditional^2 * y1)+y2))
  fscore_beta[is.na(fscore_beta)] <- 0
  fscore_beta_1 <- fscore_beta#/max(fscore_beta)
  
  df.rate_density["fscore_beta_05_normalized"] <- fscore_beta_05/max(fscore_beta_05)
  df.rate_density["fscore_beta_1_normalized"] <- fscore_beta_1/max(fscore_beta_1)
  
  df.rate_density["fscore_beta_05"] <- fscore_beta_05
  df.rate_density["fscore_beta_01"] <- fscore_beta_1
  
  ####
  
  data <- data.frame(Predictions = x,
                     Precision = y1, 
                     Recall = y2, 
                     fscore_05 = df.rate_density$fscore_beta_05_normalized,
                     fscore_1 = df.rate_density$fscore_beta_1_normalized)
  
  idx.max_05 <- max(which(data$fscore_05 == max(data$fscore_05)))
  idx.max_1 <- max(which(data$fscore_1 == max(data$fscore_1)))
  
  
  data <- data.frame(x = c(x,x,x,x),
                     y = c(y1, y2, df.rate_density$fscore_beta_05_normalized, df.rate_density$fscore_beta_1_normalized),
                     z = c(rep("Precision (number of known regulatory links)", length(x)),
                           rep("Recall (number of coregulated gene pairs with co-functional annotation)", length(x)),
                           rep("F-score (\u03B2 = 0.5)", length(x)),
                           rep("F-score (\u03B2 = 1.0)", length(x))))
  # 
  #                          rep(expression("f"[0.5]*"-score"), length(x)),
  #                          rep(expression("f"[1]*"-score"), length(x))))
  
  data$z <- factor(data$z, levels=unique(as.character(data$z)) )
  
  #cairo_pdf("example.pdf", width = 7.5, height = 6, family= "Helvetica")#"DejaVu Sans")
  
  (plot.PR <- ggplot(data, aes(x=x, y=y)) +
    geom_line(mapping=aes(colour=z,  linetype = z)) + 
    theme_bw() + 
    
    theme(legend.title=element_blank()) + 
    theme(legend.key = element_blank()) + 
    
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = c("red","steelblue","darkgray","black")) +
    scale_linetype_manual(values=c("solid", "solid","dashed", "dashed")) + 
    
    ylab("Normalized scores") +
    xlab("Fraction of predicted links / coregulated gene pairs") +  
    #scale_y_log10() +
    theme(legend.text= element_text(size=12)) + 
    theme(axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=12)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0)))
  
  #dev.off()
  
  return(list(df.rate_density=df.rate_density , plot.PR=plot.PR, ratio.gs.grn_vs_bp.normalization=ratio.gs.grn_vs_bp.normalization, max.coreg.bp=max.coreg.bp))
  
}


plot_GRACE_validation <- function(lst.eval_grace=lst.eval_grace, v.comparison_best_fscore = v.comparison_best_fscore, v.comparison_equal_size = v.comparison_equal_size, v.background = v.background){

  v.validation.grace <- unlist(lst.eval_grace)
  
  v.comparison_best_fscore <- v.comparison_best_fscore[,!names(v.comparison_best_fscore) %in% c("th.grn", "cut")]
  v.comparison_equal_size <- v.comparison_equal_size[,!names(v.comparison_equal_size) %in% c("th.grn", "cut")]
  v.background <- v.background[,!names(v.background) %in% c("th.grn", "cut")]
  
  
  n.grace <- lst.eval_grace$n.pred
  n.fscore <- v.comparison_best_fscore$n.pred
  
  # compensation vector for double annotations
  for(i in 1:length(v.validation.grace)){
    v.validation.grace[i] <- v.validation.grace[i] # - v.compensation.agris[1]  
    #v.validation.node_potentials.equal_size
    # v.validation.node_potentials.best_fscore
  }
  
  
  v.validation.grace <- c(lst.eval_grace$n.gs.grn,lst.eval_grace$n.gs.grn.Chip,lst.eval_grace$n.coreg.bp_tg_pairs,lst.eval_grace$n.coreg.localization,lst.eval_grace$n.coreg.HiC)
  v.validation.node_potentials.equal_size <- c(v.comparison_equal_size$n.gs.grn,v.comparison_equal_size$n.gs.grn.Chip,v.comparison_equal_size$n.coreg.bp_tg_pairs,v.comparison_equal_size$n.coreg.localization,v.comparison_equal_size$n.coreg.HiC)
  v.validation.node_potentials.best_fscore <- c(v.comparison_best_fscore$n.gs.grn,v.comparison_best_fscore$n.gs.grn.Chip,v.comparison_best_fscore$n.coreg.bp_tg_pairs,v.comparison_best_fscore$n.coreg.localization,v.comparison_best_fscore$n.coreg.HiC)
  

  
  v.validation.node_potentals <- c(v.background$n.gs.grn,v.background$n.gs.grn.Chip,v.background$n.coreg.bp_tg_pairs,v.background$n.coreg.localization,v.background$n.coreg.HiC)
  
  
  v.validation_background.grace <- c(lst.eval_grace$n.pred,lst.eval_grace$n.pred,lst.eval_grace$n.coreg.tg_pairs,lst.eval_grace$n.coreg.tg_pairs,lst.eval_grace$n.coreg.tg_pairs)
  names(v.validation_background.grace) <- c("regulatory_evidence", "regulatory_evidence", "cofunctional_evidence", "cofunctional_evidence", "cofunctional_evidence")
  v.validation_background.node_potentials.equal_size <- c(v.comparison_equal_size$n.pred,v.comparison_equal_size$n.pred,v.comparison_equal_size$n.coreg.tg_pairs,v.comparison_equal_size$n.coreg.tg_pairs,v.comparison_equal_size$n.coreg.tg_pairs)
  v.validation_background.node_potentials.best_fscore <- c(v.comparison_best_fscore$n.pred,v.comparison_best_fscore$n.pred,v.comparison_best_fscore$n.coreg.tg_pairs,v.comparison_best_fscore$n.coreg.tg_pairs,v.comparison_best_fscore$n.coreg.tg_pairs)
  v.validation_background.node_potentials <- c(v.background$n.pred,v.background$n.pred,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs)
  

  # customize
  validationSet = c("Gene regulatory links\n (REDfly)", 
                    "Gene regulatory links\n (ChIP binding)", 
                    "Co-functional pairs\n (Gene Ontology)",
                    "Tissue specific \n co-localized pairs\n (ImaGO)", 
                    "3D chromatin conformation\n contact pairs\n (Hi-C)")
  

#   validationSet = c("Gene regulatory links\n (AGRIS)", 
#                     "Gene regulatory links\n (ATRM)", 
#                     "Co-functional pairs\n (Gene Ontology)",
#                     "Co-localized pairs\n (Suba)", 
#                     "Metabolic pathway pairs\n (Aracyc)")
#   
  
  v.groups = c(paste("GRACE's link prioritizations (",n.grace," links)", sep = ""),
               paste("High confidence network 1 (GENIE3 + conserved motifs, ",n.grace," links)", sep = ""),
               paste("High confidence network 2 (GENIE3 + conserved motifs, ",n.fscore," links)", sep = ""))
  
  lst.res.groups <- vector(mode = "list", length = 3)
  lst.res.groups.backgrounds <- vector(mode = "list", length = 3)
  
  lst.res.groups[[1]] <- v.validation.grace
  lst.res.groups[[2]] <- v.validation.node_potentials.equal_size
  lst.res.groups[[3]] <- v.validation.node_potentials.best_fscore
  
  lst.res.groups.backgrounds[[1]] <- as.numeric(v.validation_background.grace)
  lst.res.groups.backgrounds[[2]] <- v.validation_background.node_potentials.equal_size
  lst.res.groups.backgrounds[[3]] <- v.validation_background.node_potentials.best_fscore
  
  
  #   v.groups = c(paste("GRACE (", round(lst.eval_grace$n.pred), " links)", sep = ""),
  #                paste("conserved motifs + GENIE3 (", round(lst.eval_grace$n.pred), " links)", sep = ""),
  #                paste("conserved motifs + GENIE3 (", round(lst.eval_grace$n.pred), " links)", sep = ""))
  
  
  
  df.validation <- data.frame(group = character(),
                              validationSet = character(),
                              foldChange = numeric(),
                              p.value = numeric())
  
  for(k in length(v.groups):1){
    for(j in 1:length(v.validation.node_potentals)){
      
      hitInSample <- lst.res.groups[[k]][j]
      sampleSize <- lst.res.groups.backgrounds[[k]][j]
      
      hitInPop <- v.validation.node_potentals[j]
      failInPop <- (v.validation_background.node_potentials[j] - hitInPop)
      
      foldchange <- (hitInSample/sampleSize)/(hitInPop/v.validation_background.node_potentials[j])
      
      #if(m.test == "hyper"){
        #p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
        
        
        tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
        p.val <- fisher.test(tab)$p.value
        
      #}
      
      df.validation <- rbind(df.validation, data.frame(group = v.groups[k],
                                                       validationSet = validationSet[j],
                                                       foldChange = foldchange,
                                                       p.value = p.val))
      
    }
  }
  
  df.validation["asterisk"] = ifelse(df.validation$p.val < 0.001, "**", ifelse(df.validation$p.val < 0.01, "*", ""))#, ifelse(df.validation$p.val < 0.05, "*", "")))
  

  
  #### print relative improvements #### 
  for(k in 1:length(v.validation.node_potentals)){
    
    hitInSample <- v.validation.grace[k]
    sampleSize <-  v.validation_background.grace[k]
    
    hitInPop <- v.validation.node_potentials.equal_size[k]
    popSize <- v.validation_background.node_potentials.equal_size[k]
    foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
    
    tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
    
    
    print(validationSet[k])
    print(paste("equal size comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
    
    hitInPop <- v.validation.node_potentials.best_fscore[k]
    popSize <- v.validation_background.node_potentials.best_fscore[k]
    foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
    
    tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
    fisher.test(tab)$p.value
    
    print(paste("best f score comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
    
    
    print("")
  }
  
  
    ###
    
    
    df.validation$validationSet <- factor(df.validation$validationSet, levels = df.validation$validationSet)
    
    
    data <- data.frame(Category = df.validation$group , 
                       values = df.validation$foldChange, color = c("green", "blue", "gray","red","white"), 
                       group = df.validation$validationSet)
    values_labels = as.character(df.validation$asterisk)
    
    
   # cairo_pdf("drosophila_genie3_grace.pdf", width = 7.2, height = 5, family= "Helvetica")#"DejaVu Sans")
    
      
      (barplot <- ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
        theme_bw() + 
        geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
        
        geom_text(aes(label=values_labels), fontface=2, size = 5, position=position_dodge(width=0.9),  vjust=-0.0) + 
        
        scale_fill_manual(values=c("darkgreen", "darkblue", "darkred")) + 
        
        theme(legend.title=element_blank()) + 
        
        ylab("Fold Change") +
        xlab("") +  
        
        scale_y_continuous(breaks = c(1,2,3,4,5), limits = c(-1, 5)) + 
        #scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8), limits = c(-1, 10)) + 
        
        #scale_y_log10(breaks = c(1,2,3,4,5,6,8), limits = c(-1, 10)) + 
        #labels = trans_format("log10", math_format(10^.x)))
        
        theme(legend.text=element_text(size=12)) + 
        theme(axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=10)) +
        #       theme(legend.key = element_blank()) + 
        #       
        theme(
          plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          
        ) +
        
        theme(axis.text.x=element_text(angle=45, hjust=1)) + 
        #theme( axis.text.x = element_blank()) +
        #theme(legend.background = element_rect(fill=alpha("gray90", 0.0))) +
        theme(legend.justification=c(0,1), legend.position=c(0,1)) 
      )
      
    
   # dev.off()
    
  
  
  return(barplot)
  
}


# only jaccard and th.min.coreg.tfs = 1
compute_model_assessment_parallel_optimized <- function(mat.p.grace_ensemble=mat.p.grace_ensemble, lst.benchmarks=lst.benchmarks){
            
  mat.p.crf= mat.p.grace_ensemble        
  mat.GS.grn=lst.benchmarks$mat.regulatory_evidence
  mat.bp.tgs=lst.benchmarks$mat.cofunctional_evidence
  
  
  idx.tfs <- which(rowSums(mat.p.crf) > 0)
  idx.tgs <- which(colSums(mat.p.crf) > 0)
  mat.p.crf <- mat.p.crf[idx.tfs,idx.tgs, drop = FALSE]
  
  ## atrm recovery
  tf <- intersect(rownames(mat.p.crf),rownames(mat.GS.grn)) 
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.GS.grn))
  n.gs.grn <- sum(mat.p.crf[tf,tgs] * mat.GS.grn[tf, tgs])

  #mat.grn.MR <- jaccard(t((mat.p.crf)))
  mat.grn.MR <- jaccard(t(as.matrix(mat.p.crf)))
  mat.grn.MR <- as.matrix(mat.grn.MR); 
  diag(mat.grn.MR) <- 0;
  rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(mat.p.crf) 
  mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
  
  mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
  mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
  
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.bp.tgs))
  
  mat.grn.MR <- mat.grn.MR[tgs,tgs, drop = FALSE]
  mat.bp.tgs <- mat.bp.tgs[tgs,tgs, drop = FALSE]
  
  n.coreg.tg_pairs <- sum(mat.grn.MR) 
  
  # bp pairs (subset)
  idx <- which(mat.grn.MR == 1)
  n.coreg.bp_tg_pairs <- sum(mat.bp.tgs[idx])
  
  idx <- which(colSums(mat.p.crf) > 0)
  v.tgs.set <- colnames(mat.p.crf)[idx]
  n.tg_pairs.bp <- sum(mat.bp.tgs[v.tgs.set, v.tgs.set])
  
  n.tg_pairs.bp.max <- sum(mat.bp.tgs)
  
  n.pred <- sum(mat.p.crf)
  
  rm(mat.p.crf)
  #rm(mat.grn.MR)
  #rm(mat.grn.MR.count)
  rm(mat.GS.grn)
  rm(mat.bp.tgs)
  
  return(list(n.pred=n.pred, n.gs.grn=n.gs.grn, n.coreg.bp_tg_pairs=n.coreg.bp_tg_pairs,n.coreg.tg_pairs=n.coreg.tg_pairs,
              n.tg_pairs.bp=n.tg_pairs.bp, n.tg_pairs.bp.max=n.tg_pairs.bp.max, mat.grn.MR=mat.grn.MR))
  
}



compute_grace_complete_assessment <- function(mat.grace=mat.grace, 
                                               lst.benchmarks=lst.benchmarks,
                                               lst.validation=NULL,
                                               code.validation = c("")){
                                                 
  lst.res <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.grace, lst.benchmarks=lst.benchmarks)
  
  
  lst.results <- (list(  n.pred=lst.res$n.pred,
                         n.gs.grn=lst.res$n.gs.grn,
                         
                         n.coreg.bp_tg_pairs=lst.res$n.coreg.bp_tg_pairs,
                         n.coreg.tg_pairs=lst.res$n.coreg.tg_pairs,
                         
                         n.tg_pairs.bp=lst.res$n.tg_pairs.bp,
                         n.tg_pairs.bp.max=lst.res$n.tg_pairs.bp.max))

  # subset into regulatory and co-functional evidence 
  
  if(!is.null(lst.validation)){
    
    lst.validation.results <- numeric(length(lst.validation))
    names(lst.validation.results) <- names(lst.validation)
    
    tgs <- intersect(colnames(mat.grace), colnames(lst.res$mat.grn.MR))
    mat.grn.MR <- lst.res$mat.grn.MR[tgs,tgs]
    idx <- which(mat.grn.MR == 1)
    
    tfs <- rownames(mat.grace)
    
    for(v in 1:length(lst.validation)){
      
      if(code.validation[v] == "regulatory_evidence"){
        
        lst.validation.results[v] <- sum(as.numeric(lst.validation[[v]][tfs, tgs] * mat.grace[tfs, tgs]))
        
      }else{
        
        mat.MR.validation.tgs <- lst.validation[[v]][rownames(mat.grn.MR), colnames(mat.grn.MR)]
        lst.validation.results[v] <- sum(mat.MR.validation.tgs[idx])
        
        rm(mat.MR.validation.tgs)
      }
      
    }
    
    lst.results <- c(lst.results, lst.validation.results)
  }
  
  rm(lst.res)
  
  return(lst.results)
}




compute_fmeasures_regulatory_network <- function(mat.grn=mat.gene_regulatory_network, 
                                          
                                          lst.benchmarks=lst.benchmarks,
                                          lst.validation=NULL,
                                          code.validation = c(""),
                                          
                                          n.samples = 1000,
                                          n.cpus = 2){
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.samples, n.cpus))
  registerDoParallel(cl)
  
  ### running multiple bootstrap samples in parallel
  lst.res <- foreach(i = (n.samples - 1):0, .packages=c("Matrix", "reshape2")) %dopar% {    
    
    #source("GRACE_Athaliana_evaluation.R")
    source("GRACE_helperfunctions.R")
    
    cut <- 1/n.samples * i
    cat("Processing... ", 1 - round(cut, digits = 2) , "%", "\r"); flush.console()    
    
    v.grn <- as.numeric(mat.grn[mat.grn>0])
    th.grn <- as.numeric(quantile(v.grn, cut))
    
    mat.grn.tmp <- as.matrix(mat.grn)
    mat.grn.tmp[mat.grn.tmp < th.grn] <- 0
    mat.grn.tmp[mat.grn.tmp >= th.grn] <- 1
    

    lst.results <- compute_grace_complete_assessment(mat.grace=mat.grn.tmp, 
                                                    lst.benchmarks=lst.benchmarks,
                                                    lst.validation=lst.validation,
                                                    code.validation = code.validation)
      
    lst.results <- c(lst.results, th.grn)
    lst.results <- c(lst.results, cut)
    n.idx <- length(lst.results) - 1
    names(lst.results)[n.idx:(n.idx+1)] <- c("th.grn", "cut")
    
    rm(mat.grn.tmp)
    
    return(lst.results)
  }
  stopCluster(cl)
  print(Sys.time()-strt)
  
  df.rate_density <- as.data.frame(matrix(0, nrow=n.samples, ncol = length(lst.res[[1]])))
  names(df.rate_density) <- names(lst.res[[1]])
  
  for(i in n.samples:1){
    for(j in 1:ncol(df.rate_density)){
      df.rate_density[i,j] <- lst.res[[i]][[j]]
    }
  }
    
  df.rate_density <- df.rate_density[order(df.rate_density$n.pred),]
  df.rate_density <- subset(df.rate_density, df.rate_density$th.grn != 0)
  
  return(df.rate_density)
}




### GENE REGULATORY NETWORK INFERENCE FUNCTIONS ###
infer_generegulatorynetwork <- function(mat.expression=mat.expression, method = "RAPID", return_format = "matrix", set.regulators=NULL, set.genes = NULL, seed=1234, k="sqrt", nb.trees=1000, importance.measure = "impurity", n.cpus = 2){
  
  if(method == "RAPID"){
    mat.grn <- compute_randomforest_based_GRN(mat.expression=mat.expression, k=k, nb.trees=nb.trees, set.regulators=set.regulators, set.genes = set.genes, 
                                              importance.measure=importance.measure, seed=seed, n.cpus = n.cpus)
  }
  
  if(return_format == "matrix"){
    return(mat.grn)
  }else if(return_format == "list"){
    df.grn <- as.data.frame(as.table(mat.grn), stringsAsFactors = FALSE) 
    names(df.grn)[1:3] <- c("TF","TG","v.grn")
    df.grn <- subset(df.grn, df.grn$v.grn > 0)
    df.grn <- df.grn[order(-df.grn$v.grn),]
    return(df.grn)
  }
  
}

compute_randomforest_based_GRN <- function(mat.expression=mat.expression, k="sqrt", nb.trees=1000, set.regulators = NULL, set.genes = NULL, seed=1234, importance.measure = "impurity", n.cpus = 2){
  
  library(ranger)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mat.expression.norm <- t(mat.expression)
  mat.expression.norm <- apply(mat.expression.norm, 2, function(x) { (x - mean(x)) / sd(x) } )
  n.genes <- dim(mat.expression.norm)[2]
  genes<- colnames(mat.expression.norm)
  n.samples <- dim(mat.expression.norm)[1]
  
  if(is.null(set.genes)){
    n.genes <- n.genes
    genes <- genes
    
  }else{
    n.genes <- length(set.genes)
    genes <- set.genes
    
  }
  
  #mat.expression.norm <- mat.expression.norm[,genes]
  
  if (is.null(set.regulators)) {
    n.regulators <- n.genes
    regulators <- genes
  } else {
    n.regulators <- length(set.regulators)
    # regulators provided by names
    if(is.character(set.regulators)){
      regulators <- set.regulators
      genes.undefined <- setdiff(regulators, genes)
      if (length(genes.undefined) != 0) {
        stop(paste("Error: genes: ", paste(genes.undefined, collapse = ",")," not represented in gene expression matrix \n", sep=""))
      }
      # regulators provided by indices
    }else if (is.numeric(set.regulators)) {
      regulators <- genes[set.regulators]
    }else{
      stop("Error: invalid regulator format")
    }
  }
  # set mtry
  if (class(k) == "numeric") {
    mtry <- K
  } else if (k == "sqrt") {
    mtry <- round(sqrt(n.regulators))
  } else if (k == "all") {
    mtry <- n.regulators-1
  } else {
    stop("Error: invalid parameter k, options: \"sqrt\", \"all\", or integer")
  }
  
  print(paste("Performing random-forest regression based gene regulatory network inference (# of decision trees per gene is: ", nb.trees, ", # of regulators per decision tree node: ", mtry, sep=""))
  mat.rapid <- matrix(0.0, nrow=n.regulators, ncol=n.genes, dimnames = list(regulators, genes))
  
  strt<-Sys.time()
  for(i in 1:n.genes){
    cat("Processing... ", round(i/n.genes * 100, digits = 2) , "%", "\r"); flush.console() 
    gn.i <- genes[i]
    # remove target gene from set of regulators
    regulators.i <- setdiff(regulators, gn.i)
    x <- mat.expression.norm[,regulators.i, drop=FALSE] # avoid coercion to numeric
    y <- mat.expression.norm[,gn.i]
    df.xy <- cbind(as.data.frame(x),y)
    rf.model <- ranger(y ~ ., data = df.xy, mtry=mtry, num.trees=nb.trees, importance = "impurity", num.threads = n.cpus)
    imp_scores <- rf.model$variable.importance  
    imp_scores.names <- names(imp_scores)
    mat.rapid[imp_scores.names, gn.i] <- as.numeric(imp_scores)
  }
  print(Sys.time()-strt)
  print("..finished.")
  
  return(mat.rapid / n.samples)
}       




plot_recovery_curves <- function(lst.grn.precision_per_predictions=lst.grn.precision_per_predictions,
                                 v.th.y_axis = 1.0, v.th.x_axis=1.0, 
                                 b.legend = TRUE, v.text = 12, v.text_legend = 12, 
                                 n.cores = 15, v.pos.legend = 0, v.colors = c("black", "blue", "red")){
 
  #st.legend, st.main,
  methods <- names(lst.grn.precision_per_predictions)
  
  cl<-makeCluster(n.cores)
  registerDoParallel(cl)
  lst.df.c <- foreach(i = 1:length(lst.grn.precision_per_predictions)) %dopar% { 
    #for(i in 1:length(lst.mat.grns)){
    df.c.tmp <- data.frame(lst.grn.precision_per_predictions$E.cum_precision, lst.grn.precision_per_predictions$n.pred, rep(methods, length(lst.grn.precision_per_predictions$n.predictions)))
    names(df.c.tmp) <- c("x","y","Methods")
    #lst.df.c[[i]] <- df.c.tmp
    df.c.tmp
  } 
  stopCluster(cl)
  df.c <- rbind(df.c, do.call(rbind, lst.df.c))
  
  
  plot <- (ggplot(df.c, aes(x=x, y=y, colour= Methods , fill=Methods)) +
                   geom_line() +
                   scale_color_manual(values=v.colors) + 
                   theme_bw() + 
                   ylab("Precision") +
                   xlab("Recall") +                 
                   theme(legend.text=element_text(size=v.text_legend)) + 
                   theme(legend.title=element_blank()) + 
                   theme(legend.key = element_blank()) + 
                   
                   theme(
                     plot.background = element_blank()
                     ,panel.grid.major = element_blank()
                     ,panel.grid.minor = element_blank()
                     
                   ) +
                   
                   theme(legend.background = element_rect(fill=alpha("gray90", 0.0))) +
                   
                   guides(colour = guide_legend(override.aes = list(size=2))) + 
                   
                   theme(legend.justification=c(v.pos.legend,v.pos.legend), legend.position=c(v.pos.legend,v.pos.legend)) + 
                   theme(axis.text=element_text(size=v.text, colour = "black"), axis.title=element_text(size=v.text)) + 
                   #ggtitle(paste("", st.main)) +
                   #ggtitle(paste("Precision-Recall:", st.main)) +
                   ylim(0,v.th.y_axis) +
                   xlim(0,v.th.x_axis) +
                   theme(plot.title = element_text(lineheight=1.0, face="bold", size = v.text)))
  
  return(plot=plot)
}




compute_NetworkOverlap <- function(mat.grn.prediction, mat.grn_goldstandard, trace = TRUE, test = "hyper"){
  
  tfs <- intersect(rownames(mat.grn.prediction), rownames(mat.grn_goldstandard))
  tgs <- intersect(colnames(mat.grn.prediction), colnames(mat.grn_goldstandard))
  
  mat.overlap <- mat.grn.prediction[tfs,tgs] * mat.grn_goldstandard[tfs,tgs] 
  n.pot.interactions <- nrow(mat.overlap) * ncol(mat.overlap) 
  n.grn.pred <- sum(mat.grn.prediction[tfs,tgs])
  n.grn.gs     <- sum(mat.grn_goldstandard[tfs,tgs])
  n.overlap <- sum(mat.overlap)
  n.overlap.chance <- ( n.grn.pred * n.grn.gs ) / n.pot.interactions
  fold <- n.overlap / n.overlap.chance 
  # p.val <- 1 - phyper(n.overlap-1, n.grn.pred, n.pot.interactions - n.grn.pred, n.grn.gs) # ??
  #p.val <- phyper(n.overlap-1, n.grn.pred, n.pot.interactions - n.grn.pred, n.grn.gs, lower.tail= FALSE)
  
  hitInSample <- n.overlap
  hitInPop <- n.grn.gs
  failInPop <- (n.pot.interactions - n.grn.gs)
  sampleSize <- n.grn.pred
  if(test == "hyper"){
    p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  }else if(test == "fisher"){
    # p.val <- fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize + hitInSample), 2, 2), alternative='greater')$p.value
  }
  
  if(trace){
    print(paste("predicted interactions:", n.grn.pred,", overlap:",  sum(mat.overlap), ", fold: ", fold, ", p-value ", p.val))
  }
  return(list(p.val=p.val, FoldChange=fold, networkSize=n.grn.pred, NetworkOverlap=sum(mat.overlap)))
}


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat, diag = TRUE)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat, diag = TRUE)]<- NA
  return(cormat)
}


is.installed <- function(lib)
{
  if (!require(lib,character.only = TRUE))
  {
    install.packages(lib,dep=TRUE)
    if(!require(lib,character.only = TRUE)) stop("Package not found")
  }
}


request_libraries <- function(){
  
  is.installed("graphics")
  is.installed("Matrix")
  is.installed("ggplot2")
  is.installed("ROCR")
  is.installed("reshape2")
  is.installed("CRF")
  is.installed("caTools")
  is.installed("plyr")
  is.installed("foreach")
  is.installed("doParallel")
  is.installed("igraph")
  is.installed("GenSA")
  is.installed("ranger")
  
  #detach("package:TFBSTools", unload=TRUE)
  
  if(!require(netbenchmark)){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("netbenchmark")
  }
  
  library(netbenchmark)
  library(graphics)
  library(Matrix)
  library(ggplot2)
  library(ROCR)
  library(reshape2)
  library(CRF)
  library(caTools)
  library(plyr)
  library(foreach)
  library(doParallel)
  library(igraph)
  library(GenSA)
  library(ranger)
  
}


mapping_intervals <- function(input, input_min, input_max, output_min, output_max){
  output = output_min + ((output_max - output_min) / (input_max - input_min)) * (input - input_min)
  return(output)
}



jaccard_count <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  ## only non-zero values of common
  Aim = A[im]
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim, # / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  ) 
  return( J )
}




jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  ## only non-zero values of common
  Aim = A[im]
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  ) 
  return( J )
}

