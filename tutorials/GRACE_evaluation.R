


#### BAR plots - GRACE improvements 

GRACE_improvement_analysis <- function(){

  v.methods <- names(l.methods)
  #v.methods <- c("GENIE3", "CLR","GGM", "wGLASSO", "iRafNet")
  v.groups <- c(v.methods, "GRACE")
  
  #   v.groups = c(paste("GRACE vs GRACE", sep = ""),
  #                paste("GRACE vs GENIE3", sep = ""),
  #                paste("GRACE vs CLR", sep = ""),
  #                paste("GRACE vs GGM", sep = ""),
  #                paste("GRACE vs BC3NET", sep = ""),
  #                paste("GRACE vs wGLASSO", sep = ""))
  #                # paste("GRACE vs iRafNet", sep = ""))
  #   
  # customize
  validationSet = c("Gene regulatory links\n (REDfly)", 
                    "Gene regulatory links\n (ChIP binding)", 
                    "Co-functional pairs\n (Gene Ontology)",
                    "Tissue specific \n co-localized pairs\n (ImaGO)", 
                    "3D chromatin conformation\n contact pairs\n (Hi-C)")
  
  
  validationSet <- validationSet[-1]
  v.fc.grace <- unlist(lapply(l.foldchanges.testing[[m]][[1]], function(m) m[[1]]))[-1]
  v.p.grace <- unlist(lapply(l.pvalues.testing[[m]][[1]], function(m) m[[1]]))[-1]
  
    
  # arabidopsis
#   v.fc.grace <- unlist(lapply(l.foldchanges.testing[[m]][[1]], mean))[-1]
#   v.p.grace <- unlist(lapply(l.pvalues.testing[[m]][[1]], mean))[-1]
#   
# drosophila
  
  validationSet <- validationSet[-4]
  v.fc.grace <- unlist(lapply(l.foldchanges.testing[[m]][[1]], function(m) m[[1]]))[-4]
  v.p.grace <- unlist(lapply(l.pvalues.testing[[m]][[1]], function(m) m[[1]]))[-4]
  
  
#   validationSet <- validationSet[-4]
#   v.fc.grace <- unlist(lapply(l.foldchanges.testing[[m]][[1]], mean))[-4]
#   v.p.grace <- unlist(lapply(l.pvalues.testing[[m]][[1]], mean))[-4]
#   
  
  df.validation <- data.frame(group = character(),
                              validationSet = character(),
                              foldChange = numeric(),
                              p.value = numeric())
  
  for(m in 1:length(v.methods)){
  
    v.fc.m <- unlist(lapply(l.foldchanges.testing[[m]][[2]], function(m) m[[1]]))[-4]
    v.p.m <- unlist(lapply(l.pvalues.testing[[m]][[2]], function(m) m[[1]]))[-4]
    
    #v.fc.m <- unlist(lapply(l.foldchanges.testing[[m]][[2]], mean))[-1]
    #v.p.m <- unlist(lapply(l.pvalues.testing[[m]][[2]], mean))[-1]
    #     v.fc.m <- unlist(lapply(l.foldchanges.testing[[m]][[2]], mean))[-4]
    #     v.p.m <- unlist(lapply(l.pvalues.testing[[m]][[2]], mean))[-4]
    
    for(k in 1:length(v.fc.m)){
    
      df.validation <- rbind(df.validation, data.frame(group = v.groups[m],
                                                       validationSet = validationSet[k],
                                                       foldChange = v.fc.m[k],
                                                       p.value = v.p.m[k]))
    }
  } 

  for(k in 1:length(v.fc.grace)){
    
    df.validation <- rbind(df.validation, data.frame(group = v.groups[length(v.groups)],
                                                     validationSet = validationSet[k],
                                                     foldChange = v.fc.grace[k],
                                                     p.value = v.p.grace[k]))
  }
  
  df.validation["asterisk"] = ifelse(df.validation$p.val < 0.01, "**", ifelse(df.validation$p.val < 0.05, "*", ""))#, ifelse(df.validation$p.val < 0.05, "*", "")))

  df.validation$validationSet <- factor(df.validation$validationSet, levels = df.validation$validationSet)
  
  ###
  
  for(m in 1:length(l.methods)){
    
    print(names(l.methods)[m])
    
    #v.p.grace <- unlist(lapply(l.pvalues.testing[[m]][[1]], function(m) m[[1]]))[-1]
    #l.foldchanges_comparative.testing[[m]][[2]] <- lapply(l.foldchanges_comparative.testing[[m]][[2]], function(m) m[m != "Inf"])
    
    print(lapply(l.foldchanges_comparative.testing[[m]][[2]], function(m) m[[1]]))
    print(lapply(l.pvalues_comparative.testing[[m]][[2]], function(m) m[[1]]))
    
  }
  l.foldchanges_comparative.testing <- readRDS("l.foldchanges_comparative.testing_drosophila.rds")
  l.pvalues_comparative.testing <- readRDS("l.pvalues_comparative.testing_drosophila.rds")
#   
#   
  l.foldchanges_comparative.testing <- readRDS("l.foldchanges_comparative.testing_arabidopsis.rds")
  l.pvalues_comparative.testing <- readRDS("l.pvalues_comparative.testing_arabidopsis.rds")
#   
  
  ####
  
  
  df.validation <- readRDS("df.validation_arabidopsis.rds")
  df.validation <- df.validation[c(3,4,7,8,11,12,15,16,19,20,23,24),]
  
  
  df.validation <- readRDS("df.validation_drosophila.rds")
  df.validation <- df.validation[c(2,4,6,8,10,12,14,16,18,20,22,24),]
  
  
  data <- data.frame(Category = df.validation$group , 
                     values = df.validation$foldChange, color = c("green", "blue"), #, "gray","red","white"), 
                     group = df.validation$validationSet)
  values_labels = as.character(df.validation$asterisk)
  
  #saveRDS(df.validation, "df.validation_arabidopsis.rds")
  #saveRDS(df.validation, "df.validation_drosophila.rds")
  
  #data <- data[-c(1,5,9,13,17),]
  #values_labels = values_labels[-c(1,5,9,13,17)]
  
  cairo_pdf("ensemble_validation_drosophila.pdf", width = 5, height = 5, family= "Helvetica")#"DejaVu Sans") #5
  #cairo_pdf("validation_drosophila", width = 10, height = 5, family= "Helvetica")#"DejaVu Sans") #5
  
  (ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
    theme_bw() + 
    geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
    
    geom_text(aes(label=values_labels), fontface=2, size = 5, position=position_dodge(width=0.9),  vjust=-0.0) + 
    
    scale_fill_manual(values=c("darkgray","darkorange", "darkgreen", "darkblue", "brown", "darkred")) + 
    #scale_fill_manual(values=c("darkgray","darkgreen", "darkblue", "darkred")) + 
    
    theme(legend.title=element_blank()) + 
    
    ylab("Fold Change") +
    xlab("") +  
    
    
    # arabidopsis
    # scale_y_continuous(breaks = c(1,3,5,7,9,11,13,15,17,19,21,23,25), limits = c(0, 26)) + 
    #scale_y_continuous(breaks = c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31), limits = c(0, 32)) + 
    
    # drosophila
    #scale_y_continuous(breaks = c(1,2,3,4), limits = c(0, 5)) + 
    scale_y_continuous(breaks = c(1,2,3), limits = c(0, 5)) + 
    
    #scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8), limits = c(-1, 11)) + 
    
    
    #scale_y_log10(breaks = c(1,2,3,4,5,6,8,10,15,20,25,30), limits = c(-1, 35)) + 
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
  
  dev.off()
  
  
  
}



plot_GRACE_largescale_comparative_validation <- function(lst.eval_grace=lst.eval_grace, v.comparison_best_fscore = v.comparison_best_fscore, v.comparison_equal_size = v.comparison_equal_size, v.background = v.background){
  
  ##### drosophila 
  
  v.background_datasets <- c(10, 56581, 38423, 2796802, 44335)
  
  print("characteristics of equal size comparison network")
  idx.equal_size <- max(which(df.rate_density$n.pred <= lst.eval$n.pred))
  print(df.rate_density[idx.equal_size,])
  
  df.conserved_motifs <- read.table("drosophila_datasets/motif_regnet_highestconf.txt", header  = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(df.conserved_motifs) <- c("gn.a", "gn.b", "value")
  mat.conserved_motifs <- acast(df.conserved_motifs, gn.a~gn.b, value.var = "value")
  mat.conserved_motifs[is.na(mat.conserved_motifs)] <- 0
  mat.conserved_motifs[mat.conserved_motifs != 0] <- 1
  class(mat.conserved_motifs) <- "numeric"  
  mat.conserved_motifs <- as(mat.conserved_motifs, "CsparseMatrix")

  
  lst.eval_grace=lst.eval
  
  v.validation.grace <- c(lst.eval_grace$n.gs.grn,lst.eval_grace$n.gs.grn.Chip,lst.eval_grace$n.coreg.bp_tg_pairs,lst.eval_grace$n.coreg.localization,lst.eval_grace$n.coreg.HiC)
  v.validation_background.grace <- c(lst.eval_grace$n.pred,lst.eval_grace$n.pred,lst.eval_grace$n.coreg.tg_pairs,lst.eval_grace$n.coreg.tg_pairs,lst.eval_grace$n.coreg.tg_pairs)
  names(v.validation_background.grace) <- c("regulatory_evidence", "regulatory_evidence", "cofunctional_evidence", "cofunctional_evidence", "cofunctional_evidence")

  v.background = df.rate_density[nrow(df.rate_density),]
  v.validation.node_potentals <- c(v.background$n.gs.grn,v.background$n.gs.grn.Chip,v.background$n.coreg.bp_tg_pairs,v.background$n.coreg.localization,v.background$n.coreg.HiC)
  v.validation_background.node_potentials <- c(v.background$n.pred,v.background$n.pred,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs)
  
  n.grace <- v.validation_background.grace[1]
  
  lst.res.groups <- vector(mode = "list", length = 4)
  lst.res.groups.backgrounds <- vector(mode = "list", length = 4)
  lst.res.groups[[1]] <- v.validation.grace
  lst.res.groups.backgrounds[[1]] <- as.numeric(v.validation_background.grace)
  
  
  lst.res.groups.population <- vector(mode = "list", length = 4)
  lst.res.groups.backgrounds.population <- vector(mode = "list", length = 4)
  lst.res.groups.population[[1]] <- v.validation.node_potentals
  lst.res.groups.backgrounds.population[[1]] <- as.numeric(v.validation_background.node_potentials)
  
  
  v.methods <- c("GENIE3", "CLR", "GGM")
  
  
  v.groups = c(paste("GRACE's link prioritizations", sep = ""),
               paste("GENIE3 + conserved motifs", sep = ""),
               paste("CLR + conserved motifs", sep = ""),
               paste("GGM + conserved motifs", sep = ""))
  
  for(i in 1:3){
  
    method = v.methods[i]
    
    if(FALSE){
    
      if(method == "GENIE3"){
        mat.gene_regulatory_network <- readRDS("drosophila_datasets/mat.grn.rapid.drosophila.rds")  
      }else if(method == "CLR"){
        mat.gene_regulatory_network <- readRDS("drosophila_datasets/mat.grn.clr_drosophila.rds")  
      #}else if(method == "ARACNE"){
        #mat.gene_regulatory_network <- readRDS("drosophila_datasets/mat.grn.aracne_drosophila.rds")  
      }else if (method == "GGM"){
        mat.gene_regulatory_network <- readRDS("drosophila_datasets/ggm_drosophila_0.5.rds")
      } #}else if (method == "GGM"){
        #mat.gene_regulatory_network <- readRDS("drosophila_datasets/mrnet_drosophila.rds")
      #}
      
  
      v.tgs <- intersect(colnames(mat.gene_regulatory_network), colnames(mat.conserved_motifs))
      
      v.tfs <- read.table("drosophila_datasets/tf.txt", header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
      v.tfs <- unique(v.tfs[,1])
      v.tfs <- intersect(v.tfs, rownames(mat.gene_regulatory_network))
      
      #mat.gene_regulatory_network <- mat.gene_regulatory_network[v.tfs, v.tgs]
      tfs.rapid_cns <- intersect(v.tfs, rownames(mat.conserved_motifs))
      tgs.rapid_cns <- intersect(v.tgs, colnames(mat.conserved_motifs))
      
      # 0.95 - genie3
      
      th.grn <- quantile(mat.gene_regulatory_network, 0.95) # threshold d
      # th.grn <- quantile(mat.grn.rapid, 0.95) - threshold Ath
      mat.rapid.preselection <- mat.gene_regulatory_network
      mat.rapid.preselection[mat.rapid.preselection < th.grn] <- 0
      #mat.rapid.preselection[mat.rapid.preselection > 0] <- 1
      
      mat.rapid.bound <- mat.rapid.preselection[tfs.rapid_cns, tgs.rapid_cns] * mat.conserved_motifs[tfs.rapid_cns, tgs.rapid_cns]
      mat.rapid.bound <- as.matrix(mat.rapid.bound)
      #saveRDS(mat.rapid.bound, "drosophila_datasets/mat.grn.rapid.drosophila_0.9.rds")
      
    }
    # subset 
    #   df.rate_density <- compute_fmeasures_regulatory_network(mat.grn=mat.rapid.bound, 
    #                                                           lst.benchmarks=lst.benchmarks,
    #                                                           lst.validation=lst.validation,
    #                                                           code.validation = code.validation,
    #                                                           n.samples = 1000,
    #                                                           n.cpus=n.cpus) 
  
    
    #saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_CLR_jaccard_coreg_tfs_1_1000.rds")
    #saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V2.0.rds")
    
    #saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3.0.rds")
    
    #saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V4_cell.0.rds")
    
    
    #saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_CLR_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_ARACNE_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GGM_05_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_MRNET_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    
    # load all these, subset to f0.5 score and (compute f score )
    #saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_CLR_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_ARACNE_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GGM_05_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    
    if(method == "GENIE3"){
      df.rate_density <- readRDS("drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    }else if(method == "CLR"){
      df.rate_density <- readRDS("drosophila_datasets/df.rate_density_CLR_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    }else if (method == "GGM"){
      df.rate_density <- readRDS("drosophila_datasets/df.rate_density_GGM_05_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
    }
    
    
    # fold change enrichment with respect to what?
    
    
    lst.fscore_results <- plot_fmeasure_scores(df.rate_density)
    # cairo_pdf("example.pdf", width = 7.5, height = 6, family= "Helvetica")#"DejaVu Sans")
    # lst.fscore_results$plot.PR 
    # dev.off()
    # constants needed globally
    ratio.gs.grn_vs_bp.normalization <- lst.fscore_results$ratio.gs.grn_vs_bp.normalization
    max.coreg.bp <- lst.fscore_results$max.coreg.bp
    df.rate_density <- lst.fscore_results$df.rate_density
    
    # 
    print("characteristics of equal size comparison network")
    idx.equal_size <- max(which(df.rate_density$n.pred <= lst.eval$n.pred))
    print(df.rate_density[idx.equal_size,])
    
    #     print("characteristics of best fscore comparison network")
    #     idx.best_fscore <- max(which(df.rate_density$fscore_beta_01 == max(df.rate_density$fscore_beta_01)))
    #     print(df.rate_density[idx.best_fscore,])
    
    # v.comparison_best_fscore = df.rate_density[idx.best_fscore,]
    v.comparison_equal_size = df.rate_density[idx.equal_size,]
    
    v.validation.node_potentials.equal_size <- c(v.comparison_equal_size$n.gs.grn,v.comparison_equal_size$n.gs.grn.Chip,v.comparison_equal_size$n.coreg.bp_tg_pairs,v.comparison_equal_size$n.coreg.localization,v.comparison_equal_size$n.coreg.HiC)
    #v.validation.node_potentials.best_fscore <- c(v.comparison_best_fscore$n.gs.grn,v.comparison_best_fscore$n.gs.grn.Chip,v.comparison_best_fscore$n.coreg.bp_tg_pairs,v.comparison_best_fscore$n.coreg.localization,v.comparison_best_fscore$n.coreg.HiC)
    
    v.validation_background.node_potentials.equal_size <- c(v.comparison_equal_size$n.pred,v.comparison_equal_size$n.pred,v.comparison_equal_size$n.coreg.tg_pairs,v.comparison_equal_size$n.coreg.tg_pairs,v.comparison_equal_size$n.coreg.tg_pairs)
    #v.validation_background.node_potentials.best_fscore <- c(v.comparison_best_fscore$n.pred,v.comparison_best_fscore$n.pred,v.comparison_best_fscore$n.coreg.tg_pairs,v.comparison_best_fscore$n.coreg.tg_pairs,v.comparison_best_fscore$n.coreg.tg_pairs)
  
    lst.res.groups[[1 + i]] <- v.validation.node_potentials.equal_size
    lst.res.groups.backgrounds[[1 + i]] <- as.numeric(v.validation_background.node_potentials.equal_size)
    
    
    v.background = df.rate_density[nrow(df.rate_density),]
    v.validation.node_potentals <- c(v.background$n.gs.grn,v.background$n.gs.grn.Chip,v.background$n.coreg.bp_tg_pairs,v.background$n.coreg.localization,v.background$n.coreg.HiC)
    v.validation_background.node_potentials <- c(v.background$n.pred,v.background$n.pred,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs)
    
    lst.res.groups.population[[1 + i]] <- v.validation.node_potentals
    lst.res.groups.backgrounds.population[[1 + i]] <- as.numeric(v.validation_background.node_potentials)
    
    
    
  }

  
  # customize
  validationSet = c("Gene regulatory links\n (REDfly)", 
                    "Gene regulatory links\n (ChIP binding)", 
                    "Co-functional pairs\n (Gene Ontology)",
                    "Tissue specific \n co-localized pairs\n (ImaGO)", 
                    "3D chromatin conformation\n contact pairs\n (Hi-C)")
  
  df.validation <- data.frame(group = character(),
                              validationSet = character(),
                              foldChange = numeric(),
                              p.value = numeric())
  
  for(k in length(v.groups):1){
    for(j in 1:length(v.validation.node_potentals)){
      
      hitInSample <- lst.res.groups[[k]][j]
      sampleSize <- lst.res.groups.backgrounds[[k]][j]
      
      hitInPop <- lst.res.groups.population[[k]][j]
      failInPop <- (lst.res.groups.backgrounds.population[[k]][j] - hitInPop)
      
      foldchange <- (hitInSample/sampleSize)/(hitInPop/lst.res.groups.backgrounds.population[[k]][j])
      
      tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
      p.val <- fisher.test(tab)$p.value # edit to greater instead of both
    
      df.validation <- rbind(df.validation, data.frame(group = v.groups[k],
                                                       validationSet = validationSet[j],
                                                       foldChange = foldchange,
                                                       p.value = p.val))
    }
  }
  
  df.validation["asterisk"] = ifelse(df.validation$p.val < 0.001, "**", ifelse(df.validation$p.val < 0.01, "*", ""))#, ifelse(df.validation$p.val < 0.05, "*", "")))
  
  ### 
  
  df.validation$validationSet <- factor(df.validation$validationSet, levels = df.validation$validationSet)
  
  
  data <- data.frame(Category = df.validation$group , 
                     values = df.validation$foldChange, color = c("green", "blue", "gray","red","white"), 
                     group = df.validation$validationSet)
  values_labels = as.character(df.validation$asterisk)
  
  
  cairo_pdf("drosophila_link_prioritizations.pdf", width = 7.2, height = 7, family= "Helvetica")#"DejaVu Sans") #5
  
  
  (ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
    theme_bw() + 
    geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
    
    geom_text(aes(label=values_labels), fontface=2, size = 5, position=position_dodge(width=0.9),  vjust=-0.0) + 
    
    scale_fill_manual(values=c("darkgray","darkgreen", "darkblue", "darkred")) + 
    
    theme(legend.title=element_blank()) + 
    
    ylab("Fold Change") +
    xlab("") +  
    
    scale_y_continuous(breaks = c(1,2,3,4,5,7), limits = c(-1, 7)) + 
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
  
  
  dev.off()
  
  
  

  ## all comparative heatmaps

  v.groups = c(paste("GRACE", sep = ""),
               paste("GENIE3", sep = ""),
               paste("CLR", sep = ""),
               paste("GGM", sep = ""))
  

  lst.mat.crossvalidation.foldchange <- vector(mode = "list", length = length(validationSet))
  lst.mat.crossvalidation.pvalue <- vector(mode = "list", length = length(validationSet))
  
  for(i in 1:length(validationSet)){
  
    lst.mat.crossvalidation.foldchange[[i]] <- matrix(0, nrow = length(v.groups), ncol = length(v.groups), dimnames = list(v.groups, v.groups))
    lst.mat.crossvalidation.pvalue[[i]] <- matrix(0, nrow = length(v.groups), ncol = length(v.groups), dimnames = list(v.groups, v.groups))
    
    for(k in 1:length(v.groups)){
      for(j in 1:length(v.groups)){
      
        hitInSample <- lst.res.groups[[k]][i]
        sampleSize <- lst.res.groups.backgrounds[[k]][i]
        
        hitInPop <-  lst.res.groups[[j]][i]
        failInPop <- (lst.res.groups.backgrounds[[j]][i] - hitInPop)
        
        foldchange <- (hitInSample/sampleSize)/(hitInPop/lst.res.groups.backgrounds[[j]][i])
        
        tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
        p.val <- fisher.test(tab)$p.value
        
        
        lst.mat.crossvalidation.foldchange[[i]][k,j] <- foldchange
        lst.mat.crossvalidation.pvalue[[i]][k,j] = ifelse(p.val < 0.05 & foldchange > 1, "*", "")
      
      }
    }
    
    cairo_pdf(paste("val",i,".pdf", sep = ""), width = 3.5, height = 3.5 , family= "Helvetica")#"DejaVu Sans") #5
    
    
    cormat <- lst.mat.crossvalidation.foldchange[[i]]
    melted_cormat <- melt(cormat)
    melted_cormat <- na.omit(melted_cormat)
    
    
    melted_cormat.pvalue <- melt(lst.mat.crossvalidation.pvalue[[i]])
    melted_cormat.pvalue <- na.omit(melted_cormat.pvalue)
    
    
    v.axis.ticks <- rownames(cormat)
    
    myPalette <- "black"#  ifelse(v.axis.ticks == "Non-enzyme", "gray", "black") #rainbow(length(unique(v.axis.ticks)))
    #names(myPalette) <- v.axis.ticks
    
    # Create a ggheatmap
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "black") +
      scale_fill_gradient2(low = "green", high = "red", mid = "white", 
                           midpoint = 0, limit = c(0,8), name = "fold change") +
      
      scale_x_discrete(breaks = colnames(cormat), labels=v.axis.ticks) + 
     scale_y_discrete(breaks = rownames(cormat), labels=v.axis.ticks) + 
      
      theme_minimal()+ # minimal theme
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 11, hjust = 1,
                                       colour="black")) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                       size = 11, hjust = 0,
                                       colour="black")) +
      ggtitle(validationSet[i]) +
      theme(plot.title = element_text(lineheight=1.0, face="bold", size = 16)) + 
      coord_fixed()
    
    
    
    ggheatmap <- (ggheatmap + 
                    #geom_text(aes(Var2, Var1, label = ""), color = "black", size = 4) + 
                    geom_text(aes(Var2, Var1, label = melted_cormat.pvalue$value), color = "black", size = 6) + 
                    theme(
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      panel.grid.major = element_blank()
                      #  panel.border = element_blank(),
                      #panel.background = element_blank(),
                      #axis.ticks = element_blank()
                      #legend.justification = c(1, 0),
                      #legend.position = c(0.6, 0.7),
                      #legend.direction = "horizontal")+
                      #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
                    ))
    
      plot(ggheatmap)
      
    dev.off()
      
    
    
    
    
    
    #### BAR plots - GRACE improvements 
    
    
    
    v.groups = c(paste("GRACE vs GRACE", sep = ""),
                 paste("GRACE vs GENIE3", sep = ""),
                 paste("GRACE vs CLR", sep = ""),
                 paste("GRACE vs GGM", sep = ""))
    
    
    # customize
    validationSet = c("Gene regulatory links\n (REDfly)", 
                      "Gene regulatory links\n (ChIP binding)", 
                      "Co-functional pairs\n (Gene Ontology)",
                      "Tissue specific \n co-localized pairs\n (ImaGO)", 
                      "3D chromatin conformation\n contact pairs\n (Hi-C)")
    
    df.validation <- data.frame(group = character(),
                                validationSet = character(),
                                foldChange = numeric(),
                                p.value = numeric())
    
    
    for(i in 1:length(v.validation.node_potentals)){
      for(j in 1:length(v.groups)){
        
        hitInSample <- lst.res.groups[[1]][i]
        sampleSize <- v.background_datasets[i] # lst.res.groups.backgrounds[[1]][i]
        
        hitInPop <-  lst.res.groups[[j]][i]
        failInPop <- (v.background_datasets[i] - hitInPop) #lst.res.groups.backgrounds[[j]][i]
        
        foldchange <- (hitInSample/sampleSize)/(hitInPop/v.background_datasets[i])
        #foldchange <- (hitInSample/sampleSize)/(hitInPop/lst.res.groups.backgrounds[[j]][i])
        
        tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
        p.val <- fisher.test(tab)$p.value
        
        ###
        
        df.validation <- rbind(df.validation, data.frame(group = v.groups[j],
                                                         validationSet = validationSet[i],
                                                         foldChange = foldchange,
                                                         p.value = p.val))
        
        
        
        #lst.mat.crossvalidation.foldchange[[i]][k,j] <- foldchange
        #lst.mat.crossvalidation.pvalue[[i]][k,j] = ifelse(p.val < 0.05 & foldchange > 1, "*", "")
        
      }
      
    }
    
    df.validation["asterisk"] = ifelse(df.validation$p.val < 0.01, "**", ifelse(df.validation$p.val < 0.05, "*", ""))#, ifelse(df.validation$p.val < 0.05, "*", "")))
    
    ### 
    
    df.validation$validationSet <- factor(df.validation$validationSet, levels = df.validation$validationSet)
    
    
    data <- data.frame(Category = df.validation$group , 
                       values = df.validation$foldChange, color = c("green", "blue", "gray","red","white"), 
                       group = df.validation$validationSet)
    values_labels = as.character(df.validation$asterisk)
    
    data <- data[-c(1,5,9,13,17),]
    values_labels = values_labels[-c(1,5,9,13,17)]
    
    cairo_pdf("drosophila_link_prioritizations_dataset_background.pdf", width = 7.2, height = 8, family= "Helvetica")#"DejaVu Sans") #5
    
    
    (ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
      theme_bw() + 
      geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
      
      geom_text(aes(label=values_labels), fontface=2, size = 5, position=position_dodge(width=0.9),  vjust=-0.0) + 
      
      scale_fill_manual(values=c("darkgreen", "darkblue", "darkred")) + 
      #scale_fill_manual(values=c("darkgray","darkgreen", "darkblue", "darkred")) + 
      
      theme(legend.title=element_blank()) + 
      
      ylab("Fold Change (log-scaled)") +
      xlab("") +  
      
      #scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9), limits = c(-1, 10)) + 
      #scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8), limits = c(-1, 11)) + 
      
      scale_y_log10(breaks = c(1,2,3,4,5,6,8), limits = c(-1, 12)) + 
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
      theme(legend.justification=c(1,1), legend.position=c(1,1)) 
    )
    
    
    dev.off()
    
    
    
    
    
    
    
    
      
  }
  
  
  
      
    
    
      
      
      
      
      
  }
  


evaluate_GRACE_drosophila <- function(lst.eval_grace=lst.eval_grace, v.comparison_best_fscore = v.comparison_best_fscore, v.comparison_equal_size = v.comparison_equal_size, v.background = v.background, 
                           m.test = "hyper", 
                           v.network_sizes = c("2192","792","792"), 
                           
                           v.compensation.agris = c(3,2,7), 
                           v.compensation.subacon = c(165,18,32),
                           v.compensation.aracyc = c(33,0,7)){

  mat.rapid.bound <- mat.rapid.preselection[tfs.rapid_cns, tgs.rapid_cns] * mat.cns_binding[tfs.rapid_cns, tgs.rapid_cns]
  #  saveRDS(mat.rapid.bound, "mat.rapid_de_5000.bound.rds")
  
  #####
  
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
  
  
  ### a) node potentials & relative improvement 
  
  # make 5 fold change heatmaps with values and asterisc (multiplot and gene cluster visualize)
  
  #saveRDS(lst.res.groups, "drosophila_datasets/lst.res.groups.rds")
  #saveRDS(lst.res.groups.backgrounds, "drosophila_datasets/lst.res.groups.backgrounds.rds")

  for(k in length(v.groups):1){
    for(j in 1:length(v.validation.node_potentals)){
      
      hitInSample <- lst.res.groups[[k]][j]
      sampleSize <- lst.res.groups.backgrounds[[k]][j]
      
      hitInPop <- v.validation.node_potentals[j]
      failInPop <- (v.validation_background.node_potentials[j] - hitInPop)
      
      foldchange <- (hitInSample/sampleSize)/(hitInPop/v.validation_background.node_potentials[j])
   
      tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
      p.val <- fisher.test(tab)$p.value
      
      
      df.validation <- rbind(df.validation, data.frame(group = v.groups[k],
                                                       validationSet = validationSet[j],
                                                       foldChange = foldchange,
                                                       p.value = p.val))
      
    }
  }
  
  
  
  ### c) link prioritization (compared to the maximum number of links and pairs predicted by GENIE3, the best method?)
  
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
  
  
  cairo_pdf("drosophila_genie3_grace.pdf", width = 7.2, height = 5, family= "Helvetica")#"DejaVu Sans")
  
  
  (ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
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
  
  
  dev.off()
  
  
  
  return(barplot)
  
}





#' evaluate A. thaliana enrichment GRACE vs best fscore and equal size node potential networks 
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
evaluate_GRACE <- function(lst.eval_grace=lst.eval_grace, v.comparison_best_fscore = v.comparison_best_fscore, v.comparison_equal_size = v.comparison_equal_size, v.background = v.background, 
                           m.test = "hyper", 
                           v.network_sizes = c("2192","792","792"), 
                           
                           v.compensation.agris = c(3,2,7), 
                           v.compensation.subacon = c(165,18,32),
                           v.compensation.aracyc = c(33,0,7)){
  
  ### make a dynamic version of this for publication ####
  ## multiplot - enrichment in text 
  

  v.validation.grace <- c(lst.eval_grace$n.grn_agris.gs,lst.eval_grace$n.gs.grn,lst.eval_grace$n.coreg.tg_pairs.BP,lst.eval_grace$n.coreg.tg_pairs.suba,lst.eval_grace$n.coreg.tg_pairs.aracyc)
  v.validation.grace[1] <- v.validation.grace[1] - v.compensation.agris[1]
  v.validation.grace[4] <- v.validation.grace[4] - v.compensation.subacon[1]
  v.validation.grace[5] <- v.validation.grace[5] - v.compensation.aracyc[1]
  
  v.validation_background.grace <- c(lst.eval_grace$n.pred,lst.eval_grace$n.pred,lst.eval_grace$n.coreg.tg_pairs,lst.eval_grace$n.coreg.tg_pairs,lst.eval_grace$n.coreg.tg_pairs)
  
  
  v.validation.node_potentials.equal_size <- c(v.comparison_equal_size$n.grn_agris.gs,v.comparison_equal_size$n.gs.grn,v.comparison_equal_size$n.coreg.tg_pairs.BP,v.comparison_equal_size$n.coreg.tg_pairs.suba,v.comparison_equal_size$n.coreg.tg_pairs.aracyc)
  v.validation.node_potentials.equal_size[1] <- v.validation.node_potentials.equal_size[1] - v.compensation.agris[2]
  v.validation.node_potentials.equal_size[4] <- v.validation.node_potentials.equal_size[4] - v.compensation.subacon[2]
  v.validation.node_potentials.equal_size[5] <- v.validation.node_potentials.equal_size[5] - v.compensation.aracyc[2]
  
  v.validation_background.node_potentials.equal_size <- c(v.comparison_equal_size$n.pred,v.comparison_equal_size$n.pred,v.comparison_equal_size$n.coreg.tg_pairs,v.comparison_equal_size$n.coreg.tg_pairs,v.comparison_equal_size$n.coreg.tg_pairs)
  
  
  v.validation.node_potentials.best_fscore <- c(v.comparison_best_fscore$n.grn_agris.gs,v.comparison_best_fscore$n.gs.grn,v.comparison_best_fscore$n.coreg.tg_pairs.BP,v.comparison_best_fscore$n.coreg.tg_pairs.suba,v.comparison_best_fscore$n.coreg.tg_pairs.aracyc)
  v.validation.node_potentials.best_fscore[1] <- v.validation.node_potentials.best_fscore[1] - v.compensation.agris[3]
  v.validation.node_potentials.best_fscore[4] <- v.validation.node_potentials.best_fscore[4] - v.compensation.subacon[3]
  v.validation.node_potentials.best_fscore[5] <- v.validation.node_potentials.best_fscore[5] - v.compensation.aracyc[3]
  
  v.validation_background.node_potentials.best_fscore <- c(v.comparison_best_fscore$n.pred,v.comparison_best_fscore$n.pred,v.comparison_best_fscore$n.coreg.tg_pairs,v.comparison_best_fscore$n.coreg.tg_pairs,v.comparison_best_fscore$n.coreg.tg_pairs)
 
  v.validation.node_potentals <- c(v.background$n.grn_agris.gs,v.background$n.gs.grn,v.background$n.coreg.tg_pairs.BP,v.background$n.coreg.tg_pairs.suba,v.background$n.coreg.tg_pairs.aracyc)
  v.validation_background.node_potentials <- c(v.background$n.pred,v.background$n.pred,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs,v.background$n.coreg.tg_pairs)
  
  
  validationSet = c("Gene regulatory links\n (AGRIS)", 
                    "Gene regulatory links\n (ATRM)", 
                    "Co-functional pairs\n (Gene Ontology)",
                    "Co-localized pairs\n (Suba)", 
                    "Metabolic pathway pairs\n (Aracyc)")

  v.groups = c(paste("GRACE's link prioritizations (",v.network_sizes[3]," links)", sep = ""),
               paste("High confidence network 1 (GENIE3 + conserved motifs, ",v.network_sizes[2]," links)", sep = ""),
               paste("High confidence network 2 (GENIE3 + conserved motifs, ",v.network_sizes[1]," links)", sep = ""))
  
#   v.groups = c(paste("GRACE's link prioritizations ( ",v.network_sizes[3]," links)", sep = ""),
#                paste("High confidence network 1 ( ",v.network_sizes[2]," links)", sep = ""),
#                paste("High confidence network 2 ( ",v.network_sizes[1]," links)", sep = ""))
  
  lst.res.groups <- vector(mode = "list", length = 3)
  lst.res.groups.backgrounds <- vector(mode = "list", length = 3)

  lst.res.groups[[1]] <- v.validation.grace
  lst.res.groups[[2]] <- v.validation.node_potentials.equal_size
  lst.res.groups[[3]] <- v.validation.node_potentials.best_fscore
  
  lst.res.groups.backgrounds[[1]] <- v.validation_background.grace
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
      
      if(m.test == "hyper"){
        p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
      }
      
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

                                                                                                                                                                                                    
    
  
  if(FALSE){
    
    df.validation$validationSet <- factor(df.validation$validationSet, levels = df.validation$validationSet)
    
    
    data <- data.frame(Category = df.validation$validationSet, 
                       values = df.validation$foldChange, color = c("green", "blue", "gray","red","white"), 
                       group = df.validation$group)
    values_labels = as.character(df.validation$asterisk)
    
    
    
    (barplot = ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
      theme_bw() + 
      geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
      
      geom_text(aes(label=values_labels), fontface=2, size = 5, position=position_dodge(width=0.9),  vjust=-0.0) + 
      
      scale_fill_manual(values=c("darkgreen", "darkblue", "darkred","gray","white", NA)) + 
      
      theme(legend.title=element_blank()) + 
      
      ylab("Fold Change") +
      xlab("") +  
      
      scale_y_continuous(breaks = c(1,2,3,4,5,6,8), limits = c(-1, 10)) + 
      
      #scale_y_log10(breaks = c(1,2,3,4,5,6,8), limits = c(-1, 10)) + 
      #labels = trans_format("log10", math_format(10^.x)))
      
      theme(legend.text=element_text(size=10)) + 
      theme(axis.text=element_text(size=10, colour = "black"), axis.title=element_text(size=10)) +
      #       theme(legend.key = element_blank()) + 
      #       
      theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        
      ) +
      
      #theme( axis.text.x = element_blank()) +
      #theme(legend.background = element_rect(fill=alpha("gray90", 0.0))) +
      theme(legend.justification=c(0,1), legend.position=c(0,1)) 
    )
    
  
  
  }else{
  
  ###
  
    
    df.validation$validationSet <- factor(df.validation$validationSet, levels = df.validation$validationSet)
    
    
    data <- data.frame(Category = df.validation$group , 
                       values = df.validation$foldChange, color = c("green", "blue", "gray","red","white"), 
                       group = df.validation$validationSet)
    values_labels = as.character(df.validation$asterisk)
    
    
    #cairo_pdf("example.pdf", width = 7.2, height = 5, family= "Helvetica")#"DejaVu Sans")
    (barplot <- ggplot(data=data, aes(x=factor(group), y=values, fill=Category)) + 
      theme_bw() + 
      geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
      
      geom_text(aes(label=values_labels), fontface=2, size = 5, position=position_dodge(width=0.9),  vjust=-0.0) + 
      
      scale_fill_manual(values=c("darkgreen", "darkblue", "darkred")) + 
      
      theme(legend.title=element_blank()) + 
      
      ylab("Fold Change") +
      xlab("") +  
      
      scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8), limits = c(-1, 10)) + 
      
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
    dev.off()
    
  }
  
  return(barplot)

}





#' Gene ontology visualization - as in paper
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
gene_ontology_and_phenotypoc_evaluation <- function(mat.p.grace_ensemble=mat.p.grace_ensemble, mat.grn.cns.comparison=mat.grn.cns.comparison, example = ""){
  
  
  tf.ath <- read.table("athaliana_datasets/Ath_tf_ptfb.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  v.tfs.ath.complete <- gsub("\\..*", "", tf.ath$gene_model)
  

  df.GRACE.ensemble <- as.data.frame(as.table(mat.p.grace_ensemble), stringsAsFactors = FALSE) 
  names(df.GRACE.ensemble)[1:3] <- c("TF","TG","v.grn")
  df.GRACE.ensemble <- subset(df.GRACE.ensemble, df.GRACE.ensemble$v.grn > 0)
  df.GRACE.ensemble <- df.GRACE.ensemble[,1:2]

  df.grn <- as.data.frame(as.table(mat.grn.cns.comparison), stringsAsFactors = FALSE) 
  names(df.grn)[1:3] <- c("TF","TG","v.grn")
  df.grn <- subset(df.grn, df.grn$v.grn > 0)
  df.grn <- df.grn[,1:2]
  
  library(org.At.tair.db)
  library(AnnotationDbi)
  library("GO.db")
  
  keytypes(org.At.tair.db)
  
  
  # targets
  v.tgs.grace <- unique(df.GRACE.ensemble$TG)
  v.tfs.grace <- unique(df.GRACE.ensemble$TF)
  v.tgs.grn <- unique(df.grn$TG)
  v.tfs.grn <- unique(df.grn$TF)
  v.gns.grace <- unique(c(df.GRACE.ensemble$TF, df.GRACE.ensemble$TG))
  v.gns.grn <- unique(c(df.grn$TF, df.grn$TG))
  
  
  ## only IGI and IMP 
  res.grace <- select(org.At.tair.db, v.gns.grace, "GO", "TAIR")
  res.grace <- subset(res.grace, res.grace$ONTOLOGY == "BP")
  res.grace <- subset(res.grace, res.grace$EVIDENCE %in% c("IGI", "IMP"))
  res.grace["GOTERM"] <- (select(GO.db, res.grace$GO, "TERM", "GOID"))$TERM
  res.grace.exp <- unique(res.grace[,c(1,5)])
  
  v.adapt_sets <- unique(c(res.grace.exp$GOTERM[grepl("stress", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("biotic", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("response", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("pathogen", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("fungus", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("abiotic", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("drought", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("heat", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("salt", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("cold", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("stress", res.grace.exp$GOTERM)],
                           res.grace.exp$GOTERM[grepl("adaptation", res.grace.exp$GOTERM)]))
  
  res.grace.exp <- subset(res.grace.exp, !res.grace.exp$GOTERM %in% v.adapt_sets)
  
  
  
  res.grn <- select(org.At.tair.db, v.gns.grn, "GO", "TAIR")
  res.grn <- subset(res.grn, res.grn$ONTOLOGY == "BP")
  res.grn <- subset(res.grn, res.grn$EVIDENCE %in% c("IGI", "IMP"))
  res.grn["GOTERM"] <- (select(GO.db, res.grn$GO, "TERM", "GOID"))$TERM
  res.grn.exp <- unique(res.grn[,c(1,5)])
  
  res.grn.exp <- subset(res.grn.exp, !res.grn.exp$GOTERM %in% v.adapt_sets)
  
  ###
  
  df.grace.exp <- as.data.frame(table(res.grace.exp$GOTERM))
  #df.grace.exp <- subset(df.grace.exp, df.grace.exp$Freq > 1)

  df.grn.exp <- as.data.frame(table(res.grn.exp$GOTERM))
  #df.grn.exp <- subset(df.grn.exp, df.grn.exp$Freq > 1)
  
  #res.grace.exp.only <- subset(res.grace.exp.only, !res.grace.exp.only$TAIR %in% res.grace$TAIR)
  res.grace.exp.only <- subset(res.grace.exp.only, !res.grace.exp.only$TAIR %in% res.grn.exp$TAIR)
  df.grace.exp <- as.data.frame(table(res.grace.exp.only$GOTERM))
  #df.grace.exp <- subset(df.grace.exp, df.grace.exp$Freq > 1)
  print(df.grace.exp)
  
  gns.grace.select <- res.grace.exp#[which(grepl(functional.group, res.grace.exp$GOTERM) == TRUE),]
  
  # gene names 
  df.gns_names <- read.table("athaliana_datasets/Locus_Primary_Gene_Symbol_20130117.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
  
  
  example = "auxin"
  if(example == "auxin"){

    df.GRACE.ensemble.set <- subset(df.GRACE.ensemble, df.GRACE.ensemble$TF %in% c("AT3G01330", "AT2G36010", "AT4G37750"))
    
    # example 0 - jasmonate (myc2) 
#     df.GRACE.ensemble.set <- subset(df.GRACE.ensemble, df.GRACE.ensemble$TF %in% c("AT1G32640", "AT2G38470", "AT4G34410", "AT4G17500", "AT3G44260", "AT1G74930", "AT5G47220", "AT1G24260", "AT1G59640","AT5G15800",
#                                                                                    "AT3G27810", "AT5G40350", "AT1G24260", "AT3G02310", "AT1G06160", "AT2G46510")
#                                     | df.GRACE.ensemble$TG %in% c("AT1G32640", "AT2G38470", "AT4G34410", "AT4G17500", "AT3G44260", "AT1G74930", "AT5G47220", "AT1G24260", "AT1G59640","AT5G15800",
#                                                                   "AT3G27810", "AT5G40350", "AT1G24260", "AT3G02310", "AT1G06160", "AT2G46510"))

    example = "jasmonates"
    
  }else if(example == "jasmonates"){
    
    # example 1 - cell cycle (ant)
    #df.GRACE.ensemble.set <- subset(df.GRACE.ensemble, df.GRACE.ensemble$TF %in% c("AT1G32640", "AT2G38470","AT4G25470","AT4G31800","AT2G30250","AT4G17880"))
    df.GRACE.ensemble.set <- subset(df.GRACE.ensemble, df.GRACE.ensemble$TF %in% c("AT1G32640", "AT4G17880","AT3G06740","AT5G60200","AT5G60970","AT4G25470"))
    
      
  }else{ # complete network 
    df.GRACE.ensemble.set <- df.GRACE.ensemble
  }
  
  # transform the GRACE predictions into directed graph
  g <- graph.data.frame(df.GRACE.ensemble.set, directed=TRUE)
     
  gns.grace.select.set <- subset(gns.grace.select, gns.grace.select$TAIR %in% names(V(g)))
  gns.grace.select.set <- gns.grace.select.set[-c(15),]
  
  tb.gns_functionalGroup <- table(gns.grace.select.set$TAIR)
  tb.annot_functionalGroup <- table(gns.grace.select.set$GOTERM)
  colors.go <- rainbow(length(tb.annot_functionalGroup)) # c("red", "blue", "green", "brown") #, "orange", "yellow", "pink", "cyan")
  names(colors.go) <- names(tb.annot_functionalGroup)
  
  v.shapes <- character(length(V(g)$name))
  v.color <- character(length(V(g)$name))
  
  for(i in 1:length(V(g)$name)){
    
    if(example == "auxin"){
      # example auxin
      if(V(g)$name[i] %in% c("AT2G36010", "AT3G01330")){
        v.color[i] <- "red"
      }else if(V(g)$name[i] %in% c("AT4G37750")){
        v.color[i] <- "green"
      }else{
        v.color[i] <- "darkgray"
      }
    }else{
      if(V(g)$name[i] %in% c("AT1G32640")){
        v.color[i] <- "red"
      }else if(V(g)$name[i] %in% c("AT4G17880")){
        v.color[i] <- "green"
      }else if(V(g)$name[i] %in% c("AT3G06740","AT5G60200")){
      #}else if(V(g)$name[i] %in% c("AT5G24110", "T2G38470","AT4G31800","AT2G30250", "AT4G23810", "AT4G01250", "AT4G01720", "AT2G38470")){
        v.color[i] <- "steelblue"
      }else{
        v.color[i] <- "darkgray"
      }
    }
    
    
    #     if(!is.na(as.numeric(tb.gns_functionalGroup[V(g)$name[i]]))){
    #       v.color[i] <- "green"
    #     }else{
    #       v.color[i] <- "gray"
    #     }
    #     
    if(V(g)$name[i] %in% v.tfs.ath.complete){
      v.shapes[i] = "circle"  
    }else{
      v.shapes[i] = "square"  
    }
  }
  
  # edge types
  edge.types <- numeric(length(E(g)))
  edge.colors <- numeric(length(E(g)))
  for(i in 1:length(edge.types)){
    
    tf <- names(head_of(g, E(g)[[i]]))
    tg <- names(tail_of(g, E(g)[[i]]))
    
    # check if the link is also present in the comparison network
    if(mat.grn.cns.comparison[tf,tg]){
      edge.types[i] <- 3
      edge.colors[i] <- "black"
    }else{
      edge.types[i] <- 1
      edge.colors[i] <- "red"
    }
  }
  

  node.labels <- character(length(V(g)))
  for(i in 1:length(V(g))){
    idx <- which(df.gns_names$locus_name == names(V(g))[i])
    if(length(idx) > 0){
      node.labels[i] <- df.gns_names$primary_gene_symbol[idx]
    }else{
      node.labels[i] <- names(V(g))[i]
    }
  }  

  # node names
  
  
  #l <- layout.fruchterman.reingold(g) 
  
  l <- readRDS("athaliana_datasets/layout.example_auxin.rds")
  
  #saveRDS(l, "layout.example_auxin.rds")
  #saveRDS(l, "layout.example_jasmonates.rds")
  #saveRDS(l, "layout.example_jasmonates2.rds")
  
  plot(g, 
       vertex.shape=v.shapes, 
       vertex.color=v.color,
       layout = l,
       edge.color = edge.colors, 
       edge.arrow.size=0.3,  vertex.size=3, vertex.label.cex=0.5,  edge.lty = edge.types,
       vertex.label.color = "black", vertex.label=node.labels, 
       vertex.label.font = 4,
       vertex.label.family = "Helvetica",
       edge.width=0.5,
       vertex.label.dist=0.25)

}




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
compute_grace_ensemble_model_assessment <- function(mat.p.grace_ensemble=mat.p.grace_ensemble, 
                                                    lst.benchmarks=lst.benchmarks,
                                                    mat.AtRegNet=mat.AtRegNet,
                                                    mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                    mat.aracyc.tgs=mat.aracyc.tgs,
                                                    b.jaccard = TRUE,
                                                    th.min.coreg.tfs = 1, 
                                                    beta = beta){
  
  lst.res <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.grace_ensemble,
                                             lst.benchmarks=lst.benchmarks,
                                             th.min.coreg.tfs = th.min.coreg.tfs, 
                                             b.jaccard = b.jaccard,
                                             beta = beta)
  
  
  fscore_beta=lst.res$fscore_beta
  n.pred=lst.res$n.pred
  n.gs.grn=lst.res$n.gs.grn
  
  n.coreg.bp_tg_pairs=lst.res$n.coreg.bp_tg_pairs
  n.coreg.tg_pairs=lst.res$n.coreg.tg_pairs
  
  n.tg_pairs.bp=lst.res$n.tg_pairs.bp
  n.tg_pairs.bp.max=lst.res$n.tg_pairs.bp.max
  
  mat.grn.MR = lst.res$mat.grn.MR
  
  mat.GS.grn=lst.benchmarks$mat.regulatory_evidence
  mat.bp.tgs=lst.benchmarks$mat.cofunctional_evidence
  
  tf <- intersect(rownames(mat.p.grace_ensemble),rownames(mat.AtRegNet)) 
  targets <- intersect(colnames(mat.p.grace_ensemble), colnames(mat.AtRegNet))
  
  n.grn_agris.gs <- sum(as.numeric(mat.AtRegNet[tf, targets] * mat.p.grace_ensemble[tf, targets]))
  
  # subset into regulatory and co-functional evidence 
  
  tgs <- intersect(colnames(mat.p.grace_ensemble), colnames(mat.grn.MR))
  mat.grn.MR <- mat.grn.MR[tgs,tgs]
  
  mat.MR.subacon.tgs <- mat.MR.subacon.tgs[rownames(mat.grn.MR), colnames(mat.grn.MR)]
  mat.aracyc.tgs <- mat.aracyc.tgs[rownames(mat.grn.MR), colnames(mat.grn.MR)]
  
  idx <- which(mat.grn.MR == 1)
  n.coreg.tg_pairs.suba <- sum(mat.MR.subacon.tgs[idx])
  n.coreg.tg_pairs.aracyc <- sum(mat.aracyc.tgs[idx])
  
  return(list(n.grn_agris.gs=n.grn_agris.gs,
              
              n.pred=n.pred,
              n.gs.grn=n.gs.grn,
              
              n.coreg.tg_pairs.BP=n.coreg.bp_tg_pairs,
              n.coreg.tg_pairs=n.coreg.tg_pairs,
              
              n.tg_pairs.bp=n.tg_pairs.bp,
              n.tg_pairs.bp.max=n.tg_pairs.bp.max,
              
              n.coreg.tg_pairs.suba = n.coreg.tg_pairs.suba,
              n.coreg.tg_pairs.aracyc = n.coreg.tg_pairs.aracyc))
  
}



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
compute_fmeasures_arabidopsis <- function(mat.grn=mat.gene_regulatory_network, 
                                          
                                          lst.benchmarks=lst.benchmarks,
                                          
                                          mat.AtRegNet=mat.AtRegNet,
                                          mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                          mat.aracyc.tgs=mat.aracyc.tgs,
                                          
                                          b.jaccard = TRUE,
                                          
                                          th.min.coreg.tfs = 1, 
                                          n.samples = 1000,
                                          n.cpus = 2){
  
  
  mat.GS.grn=lst.benchmarks$mat.regulatory_evidence
  mat.bp.tgs=lst.benchmarks$mat.cofunctional_evidence
  
  mat.grn.preselection <- mat.grn
  
  n.grn_agris.gs  <- numeric(n.samples)
  n.pred <- numeric(n.samples)
  n.gs.grn <- numeric(n.samples)
  n.coreg.tg_pairs.BP <- numeric(n.samples)
  n.coreg.tg_pairs <- numeric(n.samples)
  n.tg_pairs.bp <- numeric(n.samples)
  n.tg_pairs.bp.max <- numeric(n.samples)
  
  n.coreg.tg_pairs.suba <- numeric(n.samples)
  n.coreg.tg_pairs.aracyc <- numeric(n.samples)
  
  ####
  
  v.th.grn <- numeric(n.samples)
  v.th.cuts <- numeric(n.samples)
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.samples, n.cpus))
  registerDoParallel(cl)
  
  ### running multiple bootstrap samples in parallel
  lst.res <- foreach(i = (n.samples - 1):0, .packages=c("Matrix", "reshape2")) %dopar% {    
    
    source("GRACE_Athaliana_evaluation.R")
    source("GRACE_helperfunctions.R")
    
    cut <- 1/n.samples * i
    cat("Processing... ",1- round(cut, digits = 2) , "%", "\r"); flush.console()    
    
    v.grn <- as.numeric(mat.grn.preselection[mat.grn.preselection>0])
    th.grn <- as.numeric(quantile(v.grn, cut))
    
    mat.grn.tmp <- mat.grn.preselection
    mat.grn.tmp[mat.grn.tmp < th.grn] <- 0
    mat.grn.tmp[mat.grn.tmp >= th.grn] <- 1
    
    mat.p.crf=mat.p.grace_ensemble=mat.grn.tmp
    
    ######
    
    lst.res <-  compute_grace_ensemble_model_assessment(mat.p.grace_ensemble=mat.p.grace_ensemble, 
                                                        
                                                        lst.benchmarks=lst.benchmarks,
                                                        
                                                        mat.AtRegNet=mat.AtRegNet,
                                                        mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                        mat.aracyc.tgs=mat.aracyc.tgs,
                                                        
                                                        b.jaccard = b.jaccard,
                                                        
                                                        th.min.coreg.tfs = th.min.coreg.tfs, 
                                                        
                                                        beta = beta)
    
    
    lst.res <- c(lst.res, th.grn)
    lst.res <- c(lst.res, cut)
    names(lst.res)[10:11] <- c("th.grn", "cut")
    lst.res
  }
  stopCluster(cl)
  print(Sys.time()-strt)
  
    
  for(i in n.samples:1){
    
    n.grn_agris.gs[i]=lst.res[[i]]$n.grn_agris.gs
    n.gs.grn[i]=lst.res[[i]]$n.gs.grn
  
    n.pred[i]=lst.res[[i]]$n.pred
    
    n.coreg.tg_pairs.BP[i]=lst.res[[i]]$n.coreg.tg_pairs.BP
    n.coreg.tg_pairs[i]=lst.res[[i]]$n.coreg.tg_pairs
    
    n.tg_pairs.bp[i]=lst.res[[i]]$n.tg_pairs.bp
    n.tg_pairs.bp.max[i]=lst.res[[i]]$n.tg_pairs.bp.max
    
    n.coreg.tg_pairs.suba[i] = lst.res[[i]]$n.coreg.tg_pairs.suba
    n.coreg.tg_pairs.aracyc[i] = lst.res[[i]]$n.coreg.tg_pairs.aracyc
    
    v.th.grn[i] <- lst.res[[i]]$th.grn
    v.th.cuts[i] <- round(lst.res[[i]]$cut, digits = 5)
  }
  
    
    
  
  # load dataset
  df.rate_density <- data.frame(n.grn_agris.gs=n.grn_agris.gs,
                             
                                n.pred=n.pred,
                                n.gs.grn=n.gs.grn,
                                n.coreg.tg_pairs.BP=n.coreg.tg_pairs.BP,
                                n.coreg.tg_pairs=n.coreg.tg_pairs,
                                
                                n.tg_pairs.bp=n.tg_pairs.bp,
                                n.tg_pairs.bp.max=n.tg_pairs.bp.max,
                                
                                n.coreg.tg_pairs.suba=n.coreg.tg_pairs.suba,
                                n.coreg.tg_pairs.aracyc=n.coreg.tg_pairs.aracyc,
                      
                                v.th.grn=v.th.grn,
                                v.th.cuts=v.th.cuts)
  
  df.rate_density <- df.rate_density[order(df.rate_density$n.pred),]
  df.rate_density <- subset(df.rate_density, df.rate_density$v.th.grn != 0)
  
  return(df.rate_density)
}




#' load example datasets for A. thaliana development network inference and evaluation
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param b.jaccard definition of coregulated gene pairs.
#' @param beta definition of coregulated gene pairs.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
compute_fscore_model_assessment <- function(mat.p.grace_ensemble=mat.p.grace_ensemble,
                                            lst.benchmarks=lst.benchmarks,
                                            th.min.coreg.tfs = 1, 
                                            b.jaccard = TRUE,
                                            beta = 1.0){
  
  mat.p.crf= mat.p.grace_ensemble        
  mat.GS.grn=lst.benchmarks$mat.regulatory_evidence
  mat.bp.tgs=lst.benchmarks$mat.cofunctional_evidence

  ## atrm recovery
  tf <- intersect(rownames(mat.p.crf),rownames(mat.GS.grn))
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.GS.grn))
  n.gs.grn <- sum(mat.p.crf[tf,tgs] * mat.GS.grn[tf, tgs])
  
  ### coop-coreg tg pairs
  
  mat.grn.MR.count <- jaccard_count(t((mat.p.crf)))
  #mat.grn.MR.count <- jaccard_count(t(as.matrix(mat.p.crf)))
  #mat.grn.MR.count <- as.matrix(mat.grn.MR.count); 
  diag(mat.grn.MR.count) <- 0;
  rownames(mat.grn.MR.count) <- colnames(mat.grn.MR.count) <- colnames(mat.p.crf) 
  mat.grn.MR.count <- as(mat.grn.MR.count, "CsparseMatrix")
  
  mat.grn.MR.count@x[mat.grn.MR.count@x < th.min.coreg.tfs] <- 0
  mat.grn.MR.count@x[mat.grn.MR.count@x >= th.min.coreg.tfs] <- 1
  
  if(b.jaccard){ # jaccard > 0.5 rule
    
    mat.grn.MR <- jaccard(t((mat.p.crf)))
    #mat.grn.MR <- jaccard(t(as.matrix(mat.p.crf)))
    #mat.grn.MR <- as.matrix(mat.grn.MR); 
    diag(mat.grn.MR) <- 0;
    rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(mat.p.crf) 
    mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
 
    mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
    mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
    
    mat.grn.MR.count <- mat.grn.MR * mat.grn.MR.count
  }
  
  
  
  
  #mat.grn.MR.count[lower.tri(mat.grn.MR.count)] <- 0 
  mat.grn.MR <- mat.grn.MR.count
  
  tgs <- intersect(colnames(mat.p.crf), colnames(mat.bp.tgs))
  
  mat.grn.MR <- mat.grn.MR[tgs,tgs]
  mat.bp.tgs <- mat.bp.tgs[tgs,tgs]
  
  n.coreg.tg_pairs <- sum(mat.grn.MR) 
  
  # bp pairs (subset)
  idx <- which(mat.grn.MR == 1)
  n.coreg.bp_tg_pairs <- sum(mat.bp.tgs[idx])
  
  idx <- which(colSums(mat.p.crf) > 0)
  v.tgs.set <- colnames(mat.p.crf)[idx]
  n.tg_pairs.bp <- sum(mat.bp.tgs[v.tgs.set, v.tgs.set])
  
  n.tg_pairs.bp.max <- sum(mat.bp.tgs)
  
  n.pred <- sum(mat.p.crf)
  
  
  return(list(n.pred=n.pred, n.gs.grn=n.gs.grn, n.coreg.bp_tg_pairs=n.coreg.bp_tg_pairs,n.coreg.tg_pairs=n.coreg.tg_pairs,
              n.tg_pairs.bp=n.tg_pairs.bp, n.tg_pairs.bp.max=n.tg_pairs.bp.max, mat.grn.MR=mat.grn.MR))
              
}









#' perform A. thaliana gene annotation
#' 
#' @param mat.grn A number.
#' @param mat.FunctionalAssociation A number.
#' @param th.preselection .
#' @param n.cores Number of cores to use.
#' @return adj.links \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
perform_Athaliana_functional_gene_annotation <- function(mat.grn.selection=mat.grn.selection,  keytype = keytype, annotation = "GO", go.type = c("MF", "BP", "CC"), evidences = c("IGI", "IMP"), min.genes_per_annotation = 3, th.p.value = 0.01, 
                                                         legend = FALSE, v.filter = c("regulation", "defense", "stress", "biotic", "pathogen", "fungus", "abiotic", "drought", "heat", "salt", "cold", "adaptation")){
  
  # online enrichment analysis of the grace network 
  
  
  is.installed("AnnotationDbi")
  is.installed("GO.db")
  
  library(AnnotationDbi)
  library("GO.db")
  
  
  if(!require(org.At.tair.db)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.At.tair.db")
  }
  
  library(org.At.tair.db)
  keytypes(org.At.tair.db)
  AnnotationDb = org.At.tair.db
  
  ## extract gene ontology information for entire top prediction network 
  idx.rows <- which(rowSums(mat.grn.selection) > 0)
  idx.cols <- which(colSums(mat.grn.selection) > 0)
  mat.grn.selection <- mat.grn.selection[idx.rows, idx.cols]
  
  mat.grn.selection[mat.grn.selection > 0] <- 1
  #gns <- unique(colnames(mat.grn.selection))
  gns <- unique(c(rownames(mat.grn.selection), colnames(mat.grn.selection)))
  
  n.genome <- length(keys(AnnotationDb))
  res.go.genome <- select(AnnotationDb, keys(AnnotationDb), annotation, keytype=keytype)
  
  res.go.genome <- subset(res.go.genome, res.go.genome$EVIDENCE %in% evidences)
  res.go.genome["GOTERM"] <- (select(GO.db, res.go.genome$GO, "TERM", "GOID"))$TERM
  #res.go <- unique(res.go[,c(1,5)])
  
  tb.go.genome <- as.data.frame(table(res.go.genome$GOTERM))
  tb.go.genome <- subset(tb.go.genome, tb.go.genome$Freq >= min.genes_per_annotation)
  res.go.genome <- subset(res.go.genome, res.go.genome$GOTERM %in% as.character(tb.go.genome[,1]))
  
  res.go.grn <- select(AnnotationDb, gns, annotation, keytype=keytype) 
  res.go.grn <- subset(res.go.grn, res.go.grn$ONTOLOGY == go.type)
  res.go.grn <- subset(res.go.grn, res.go.grn$EVIDENCE %in% evidences)
  res.go.grn["GOTERM"] <- (select(GO.db, res.go.grn$GO, "TERM", "GOID"))$TERM
  #res.go <- unique(res.go[,c(1,5)])
  
  tb.go.grn <- as.data.frame(table(res.go.grn$GOTERM))
  tb.go.grn <- subset(tb.go.grn, tb.go.grn$Freq >= min.genes_per_annotation)
  res.go.grn <- subset(res.go.grn, res.go.grn$GOTERM %in% as.character(tb.go.grn[,1]))
  
  n.gns.grn <- length(gns)
  
  tb.go.grn["foldchange"] <- 0
  tb.go.grn["p.val"] <- 1
  names(tb.go.grn) <- c("GOTERM", "Freq", "foldchange", "p.val")
  
  #go.terms <- table(res.go.c$GOTERM)
  for(i in 1:nrow(tb.go.grn)){
    
    hitInSample <- as.numeric(tb.go.grn$Freq[i])
    sampleSize <- n.gns.grn
    
    idx <- which(as.character(tb.go.grn$GOTERM[i]) == tb.go.genome$Var1)
    hitInPop <- tb.go.genome[idx,2]
    failInPop <- n.genome - hitInPop
    
    # barplots per module?  - make hypergeometric test
    #print((hitInSample/sampleSize) / (hitInPop / length(gns)))
    tab <- matrix(c(hitInSample, hitInPop, sampleSize, length(gns) - hitInPop), nrow = 2)
    
    tb.go.grn$foldchange[i] <- ((hitInSample/sampleSize) / (hitInPop / n.genome))
    tb.go.grn$p.val[i] <- p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE) #fisher.test(tab)$p.value
  }
  
  tb.go.grn <- subset(tb.go.grn, tb.go.grn$p.val <= th.p.value)
  names(tb.go.grn) <- c("GOTERM", "Freq", "foldchange", "p.val")
  
  
  ### FILTER BPs ###
  
  if(length(v.filter) > 0){
    v.adapt_sets <- character()
    
    for(j in 1:length(v.filter)){
      v.adapt_sets <- c(v.adapt_sets, unique(as.character(tb.go.grn$GOTERM[grepl(v.filter[j], tb.go.grn$GOTERM)])))
    }
    
    tb.go.grn <- subset(tb.go.grn, !tb.go.grn$GOTERM %in% v.adapt_sets)
  }
  
  ###
  
  data <- data.frame(Category = tb.go.grn$GOTERM, 
                     p.val = tb.go.grn$p.val,
                     values = tb.go.grn$foldchange)
  data <- data[nrow(data):1,]
  data$Category <- factor(data$Category, levels = data$Category)
  
  #group = tb.go.grn$set)
  
  data$p.val = ifelse(data$p.val < 0.001, "**", ifelse(data$p.val < 0.01, "*", ""))
  # values_labels = ifelse(tb.go.grn$p.val <= 0.001, "***", ifelse(tb.go.grn$p.val <= 0.01, "**", ifelse(tb.go.grn$p.val <= 0.05, "*", "")))
  
  
  if(legend){
    
    (barplot = ggplot(data=data, aes(x=factor(Category), y=values, fill=Category)) + 
       theme_bw() + 
       geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
       
       geom_text(aes(label=values_labels), position=position_dodge(width=0.9),  vjust=-0.05) + 
       
       theme(legend.title=element_blank()) + 
       theme(legend.key = element_blank()) + 
       
       theme(plot.background = element_blank()
             ,panel.grid.major = element_blank()
             ,panel.grid.minor = element_blank()) +
       
       ylab("Fold Change") +
       xlab("") +  
       scale_y_log10() +
       theme(legend.text=element_text(size=12)) + 
       
       theme(axis.text=element_text(size=12, colour = "black"), axis.title=element_text(size=12)) +
       theme( axis.text.x = element_blank()))
    
  }else{
    
    (barplot = ggplot(data=data, aes(x=factor(Category), y=values, fill=Category)) + 
      theme_bw() + 
      geom_bar(aes(fill = Category), stat="identity", position = "dodge",color="black") + 
      #geom_text(aes(label=values_labels), position=position_dodge(width=0.9),  vjust=-0.05) + 
      geom_text(aes(label= data$p.val), position=position_dodge(width=0.9),  hjust=-0.15, vjust = 0.7) + 
      theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      
      ylab("Fold Change") +
      xlab("") +  
      
      scale_y_sqrt(breaks = c(5,10,15,20,30), limits = c(0,40)) + #(breaks = c(25,50,100,150,250,350,500), limits = c(1,600)) + 
      coord_flip() +
      
      scale_fill_manual(values = rep("steelblue", nrow(data))) +
      #scale_fill_grey(start = 0, end = .9) + 
      
      theme(axis.text.x=element_text(angle=45, hjust=1)) + 
      
      theme(axis.text=element_text(size=8, colour = "black"), axis.title=element_text(size=8)) + 
      guides(fill=FALSE))
    
    
  }
  
  
  # perform a fisher exact test with respect to complete table
  names(tb.go.grn) <- c("annotation", "numberOfgenes", "foldchange", "p.value")
  res.go.grn <- subset(res.go.grn , res.go.grn$GOTERM %in% tb.go.grn$annotation)
  
  return(list(go.annotations=res.go.grn, go.frequency_table=tb.go.grn, go.enrichment_barplot=barplot))
}



