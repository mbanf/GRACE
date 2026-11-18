


print("Computing averaged bootstrap error validation - using same size network")
### compute the bootstrap averaged errors 

l.foldchanges <- vector(mode = "list", length = 5)
names(l.foldchanges) <- c("GENIE3", "CLR",  "GGM", "WGLASSO", "iRafNet")

for(m in 1:6){
  l.foldchanges[[m]] <- vector(mode = "list", length = 2)
  names(l.foldchanges[[m]]) <- c("GRACE", "Original")
  for(j in 1:2){
    l.foldchanges[[m]][[j]] <- vector(mode = "list", length = 5)
    for(s in 1:5){
      l.foldchanges[[m]][[j]][[s]] <- numeric(n.models)
    }
  }
}

l.pvalues <- vector(mode = "list", length = 5)
names(l.pvalues) <- c("GENIE3", "CLR", "GGM", "WGLASSO", "iRafNet")

for(m in 1:6){
  l.pvalues[[m]] <- vector(mode = "list", length = 2)
  names(l.pvalues[[m]]) <- c("GRACE", "Original")
  for(j in 1:2){
    l.pvalues[[m]][[j]] <- vector(mode = "list", length = 5)
    for(s in 1:5){
      l.pvalues[[m]][[j]][[s]] <- numeric(n.models)
    }
  }
}

l.foldchanges_comparative <- vector(mode = "list", length = 5)
names(l.foldchanges_comparative) <- c("GENIE3", "CLR",  "GGM", "WGLASSO", "iRafNet")

for(m in 1:6){
  l.foldchanges_comparative[[m]] <- vector(mode = "list", length = 2)
  names(l.foldchanges_comparative[[m]]) <- c("GRACE", "Original")
  for(j in 1:2){
    l.foldchanges_comparative[[m]][[j]] <- vector(mode = "list", length = 5)
    for(s in 1:5){
      l.foldchanges_comparative[[m]][[j]][[s]] <- numeric(n.models)
    }
  }
}

l.pvalues_comparative <- vector(mode = "list", length = 5)
names(l.pvalues_comparative) <- c("GENIE3", "CLR",  "GGM", "WGLASSO", "iRafNet")

for(m in 1:6){
  l.pvalues_comparative[[m]] <- vector(mode = "list", length = 2)
  names(l.pvalues_comparative[[m]]) <- c("GRACE", "Original")
  for(j in 1:2){
    l.pvalues_comparative[[m]][[j]] <- vector(mode = "list", length = 5)
    for(s in 1:5){
      l.pvalues_comparative[[m]][[j]][[s]] <- numeric(n.models)
    }
  }
}

validationSet = c("Gene regulatory links\n (AGRIS)", 
                  "Gene regulatory links\n (ATRM)", 
                  "Co-functional pairs\n (Gene Ontology)",
                  "Co-localized pairs\n (Suba)", 
                  "Metabolic pathway pairs\n (Aracyc)")

v.sets_map <- c(1,3,4,8,9)
#v.sets_map <- c(1,6,7,15,16)
names(v.sets_map) <- validationSet


validationSet = c("Gene regulatory links\n (REDfly)", 
                  "Gene regulatory links\n (ChIP binding)", 
                  "Co-functional pairs\n (Gene Ontology)",
                  "Tissue specific \n co-localized pairs\n (ImaGO)", 
                  "3D chromatin conformation\n contact pairs\n (Hi-C)")

v.sets_map <- c(2,7,3,8,9)
names(v.sets_map) <- validationSet

v.background = df.rate_density[nrow(df.rate_density),]
n.pred.background <- df.rate_density$n.pred[nrow(df.rate_density)]
n.coreg.background <- df.rate_density$n.coreg.tg_pairs[nrow(df.rate_density)]
v.backgroundSet <- as.numeric(df.rate_density[nrow(df.rate_density),v.sets_map]) # needs to be adjusted
names(v.backgroundSet) <- validationSet

# generate training and test bootstrap prediction errors

# GENIE3, GGM, CLR, BC3NET, wGLASSO, iRafNet
l.methods <- vector(mode = "list", length = 5)


m.grace_initial <- as.matrix(readRDS("drosophila_datasets/mat.grn.rapid.drosophila.rds"))
v.tfs <- rownames(m.grace_initial)
v.tgs <- colnames(m.grace_initial)


l.methods[[1]] <- readRDS("athaliana_datasets/mat.conserved_motifs_plus_genie3_development.rds")
#l.methods[[1]] <- as.matrix(readRDS("drosophila_datasets/mat.grn.rapid.drosophila.rds"))#[v.tfs,] # genie 3
l.methods[[1]] <- as.matrix(readRDS("drosophila_datasets/randomforest.rds"))# genie 3
v.tfs <- rownames(l.methods[[1]])

l.methods[[2]] <- as.matrix(readRDS("drosophila_datasets/clr.rds")[v.tfs,])

#l.methods[[2]] <- as.matrix(readRDS("drosophila_datasets/mat.grn.clr_drosophila.rds")[v.tfs,])
#l.methods[[3]] <- readRDS("drosophila_datasets/mat.grn.aracne_drosophila.rds")[v.tfs,v.tgs]
# l.methods[[4]] <- readRDS("bc3net_drosophila.rds")[v.tfs,v.tgs]
l.methods[[3]] <- as.matrix(readRDS("drosophila_datasets/ggm_drosophila_0.5.rds")[v.tfs,])
#l.methods[[5]] <- 

# glasso <- readRDS("glasso_drosophila.rds")
# rownames(glasso$wi) <- colnames(glasso$wi) <- v.gns
# saveRDS(glasso$wi[v.tfs,], "drosophila_datasets/wGlasso_drosophila.rds")

l.methods[[4]] <- as.matrix(abs(readRDS("drosophila_datasets/wGlasso_drosophila.rds")[v.tfs,]))

# l.methods[[4]] <- matrix(0, nrow = length(v.tfs), ncol = length(v.tgs), dimnames = list(v.tfs, v.tgs))
# m.tmp <- as.matrix(abs(readRDS("drosophila_datasets/wGlasso_drosophila.rds")[v.tfs,]))
# l.methods[[4]][v.tfs, intersect(v.tgs, colnames(m.tmp))] <- m.tmp[v.tfs, intersect(v.tgs, colnames(m.tmp))]


l.methods[[5]] <- as.matrix(readRDS("iRafNet_drosophila.rds")[v.tfs,])

##

l.methods[[1]] <- readRDS("athaliana_datasets/mat.rapid_de_5000.rds") #mat.conserved_motifs_plus_genie3_development.rds") # genie 3 
v.tfs <- rownames(l.methods[[1]])
v.tgs <- colnames(l.methods[[1]])

l.methods[[2]] <- readRDS("athaliana_datasets/clr.rds")[v.tfs,] # no binding
l.methods[[3]] <- readRDS("inferred.pcor_arabidopsis.rds")[v.tfs,] # no binding 
# m.tmp <- readRDS("glasso_arabidopsis.rds")
# rownames(m.tmp$wi) <- colnames(m.tmp$wi) <- v.gns
# saveRDS(m.tmp$wi[v.tfs,], "athaliana_datasets/wglasso_arabidopsis.rds")

l.methods[[4]]  <- abs(readRDS("athaliana_datasets/wglasso_arabidopsis.rds")[v.tfs,])
l.methods[[5]]  <- abs(readRDS("iRafNet_arabidopsis.rds")[v.tfs,])

for(i in 1:5){
  print(i)
  l.methods[[i]] <- l.methods[[i]][,v.tgs]
  th <- quantile(l.methods[[i]], 0.95)
  l.methods[[i]][l.methods[[i]] < th] <- 0
}
# 


for(i in 1:5){
  print(i)
  l.methods[[i]] <- l.methods[[i]][v.tfs,]
}
#saveRDS(m.clr[v.tfs,v.tgs], "athaliana_datasets/mat.clr_arabidopsis.rds")

#l.methods[[3]] <- readRDS("drosophila_datasets/mat.grn.aracne_drosophila.rds")[v.tfs,v.tgs]
# l.methods[[4]] <- readRDS("bc3net_drosophila.rds")[v.tfs,v.tgs]
#l.methods[[3]] <- readRDS("athaliana_datasets/drosophila_datasets/ggm_drosophila.rds")[v.tfs,v.tgs]



#l.methods[[5]] <- readRDS("glasso_drosophila.rds")
# l.methods[[4]] <- abs(readRDS("glasso_arabidopsis.rds")[v.tfs,v.tgs])
#l.methods[[4]] <- readRDS("glasso_arabidopsis.rds")


#m.tmp <- readRDS("glasso_drosophila.rds")
#l.methods[[4]]  <- abs(readRDS("athaliana_datasets/wglasso_arabidopsis.rds")[v.tfs,])

# rownames(l.methods[[4]]$wi) <- colnames(l.methods[[4]]$wi) <- v.gns
#saveRDS(l.methods[[4]]$wi[v.tfs,v.tgs], "athaliana_datasets/wglasso_arabidopsis.rds")

#rownames(l.methods[[5]]$wi) <- colnames(l.methods[[5]]$wi) <- v.gns
#saveRDS(l.methods[[5]]$wi[v.tfs,v.tgs], "drosophila_datasets/wGlasso_drosophila.rds")
# l.methods[[5]]
# names(l.methods) <- c("GENIE3", "CLR", "GGM", "wGLASSO", "iRafNet")
# 
# # readRDS("iRafNet_drosophila")
# 
# lst.validation <- list(lst.data$mat.AtRegNet, lst.data$mat.MR.subacon.tgs, lst.data$mat.aracyc.tgs)
# # v.gs.copairs <- c(sum(lst.validation[[1]][v.tfs,v.tgs]), sum(lst.validation[[2]][v.tgs,v.tgs]), sum(lst.validation[[3]][v.tgs,v.tgs]))
# 
# #v.gs.copairs <- c(sum(lst.validation$n.gs.grn.Chip[v.tfs,v.tgs]), sum(lst.validation$n.coreg.localization[v.tgs,v.tgs]), sum(lst.validation$n.coreg.HiC[v.tgs,v.tgs]))
# 
# ####
# 
# for(t in 1:2){
#   
#   if(t == 1){
#     b.validation <- FALSE
#   }else{
#     b.validation <- TRUE
#   }
#   
#   for(i in 1:n.models){
#     
#     cat("Processing... ", round(i/n.models * 100, digits = 2) , "%", "\r"); flush.console() 
#     
#     # use different seeds per run (simulated annealing)
#     v.seed = 1234 + 158 * i
#     set.seed(v.seed)
#     
#     #### Bootstrapping (based on regulators based meta module sampling (0.632))
#     tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
#     tgs.mrf <- unique(lst.grace_support$df.grn$TG)
#     
#     lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size, b.validation = b.validation)
#     lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
#     lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
#     
#     lst.benchmarks <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
#     names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
#     
#     lst.validation <- list(lst.data$mat.Chip, lst.data$mat.tissue_localization.tgs, lst.data$mat.HiC)
#     names(lst.validation) <- c("n.gs.grn.Chip","n.coreg.localization","n.coreg.HiC")
#     code.validation <- c("regulatory_evidence", "cofunctional_evidence", "cofunctional_evidence")
#     
#     ###
#     
#     alpha  <- lst.grace_models[[i]]$weights[1]
#     gamma  <- lst.grace_models[[i]]$weights[2]
#     lambda <- lst.grace_models[[i]]$weights[3]
#     lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields(lst.grace_support=lst.grace_support, lst.modules=lst.modules, 
#                                                               v.alpha = alpha, v.gamma=gamma, v.lambda=lambda, b.prim=TRUE) 
#     mat.p.crf <- lst.mat.p.crf$mat.p.crf
#     
#     ### 
#     
#     l.eval <- vector(mode = "list", length = 2)
#     
#     l.eval[[1]] <- compute_grace_complete_assessment(mat.grace=mat.p.crf, 
#                                                      lst.benchmarks=lst.benchmarks,
#                                                      lst.validation=lst.validation,
#                                                      code.validation = code.validation)
#     
#     
#     n.links <- sum(mat.p.crf)
#     
#     
#     # background links and copais - similar for all 
#     n.links_background <- 1414 * 22591 # global  #length(v.tfs) * length(v.tgs)
#     n.co_pairs <- 22591 * 22591 #length(v.tgs) * length(v.tgs)
#     
#     
#     n.gs_background <- sum(mat.regulatory_evidence)
#     
#     
#     
#     
#     v.gs.grn <- c(n.gs_background, n.chip_background)
#     
#     # all methods as comparison models here #
#     for(m in 1:length(l.methods)){
#       
#       mat.gene_regulatory_network <- as.matrix(l.methods[[m]])
#       
#       mat.grn.comparison <- mat.gene_regulatory_network[v.tfs,v.tgs]
#       mat.grn.comparison <- mat.grn.comparison[rownames(lst.benchmarks$mat.regulatory_evidence),]
#       
#       df.grn.comparison <- as.data.frame(as.table(mat.grn.comparison), stringsAsFactors = FALSE) 
#       names(df.grn.comparison)[1:3] <- c("TF","TG","v.grn")
#       
#       
#       df.grn.comparison <- subset(df.grn.comparison, df.grn.comparison$v.grn > 0)
#       
#       
#       df.grn.comparison <- df.grn.comparison[order(-df.grn.comparison$v.grn),] 
#       df.grn.comparison <- df.grn.comparison[1:n.links,]
#       
#       mat.grn.comparison <- acast(df.grn.comparison, TF~TG, value.var="v.grn")
#       mat.grn.comparison[is.na(mat.grn.comparison)] <- 0; 
#       mat.grn.comparison[mat.grn.comparison > 0] <- 1
#       
#       l.eval[[2]] <- compute_grace_complete_assessment(mat.grace=mat.grn.comparison, 
#                                                        lst.benchmarks=lst.benchmarks,
#                                                        lst.validation=lst.validation,
#                                                        code.validation = code.validation)
#       
#       for(j in 1:2){
#         
#         # print(j)
#         for(s in 1:length(v.sets_map)){
#           
#           hitInSample <- l.eval[[j]][[v.sets_map[s]]]
#           
#           if(s <= 2){
#             hitInPop <- v.gs.grn[s]
#             sampleSize <- n.links
#             failInPop <- (n.links_background - hitInPop)
#             foldchange <- (hitInSample/sampleSize)/(hitInPop/n.links_background)
#           }else{
#             hitInPop <- v.gs.copairs[s - 2]
#             sampleSize <- l.eval[[j]]$n.coreg.tg_pairs
#             failInPop <- (n.co_pairs - hitInPop)
#             foldchange <- (hitInSample/sampleSize)/(hitInPop/n.co_pairs)
#           }
#           
#           #print(foldchange)
#           
#           p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
#           
#           l.foldchanges[[m]][[j]][[s]][[i]] <- foldchange
#           l.pvalues[[m]][[j]][[s]][[i]] <- p.val
#           
#         }
#         #print("")
#       }
#       
#       ### comparative enrichment 
#       
#       for(s in 1:length(v.sets_map)){
#         
#         hitInSample <- l.eval[[1]][[v.sets_map[s]]]
#         hitInPop <- l.eval[[2]][[v.sets_map[s]]]
#         
#         if(s <= 2){
#           sampleSize <- n.links
#           failInPop <- (n.links - hitInPop)
#           foldchange <- (hitInSample/n.links)/(hitInPop/n.links)
#         }else{
#           sampleSize <- l.eval[[1]]$n.coreg.tg_pairs
#           failInPop <- (l.eval[[2]]$n.coreg.tg_pairs - hitInPop)
#           foldchange <- (hitInSample/sampleSize)/(hitInPop/l.eval[[2]]$n.coreg.tg_pairs)
#         }
#         
#         tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
#         p.val <- fisher.test(tab)$p.value
#         
#         l.foldchanges_comparative[[m]][[j]][[s]][[i]] <- foldchange
#         l.pvalues_comparative[[m]][[j]][[s]][[i]] <- p.val
#         
#       }
#     }
#   }
#   
#   if(t == 1){
#     l.foldchanges.training <- l.foldchanges
#     l.pvalues.training <- l.pvalues
#     l.foldchanges_comparative.training <- l.foldchanges_comparative
#     l.pvalues_comparative.training <- l.pvalues_comparative
#   }else{
#     l.foldchanges.testing <- l.foldchanges
#     l.pvalues.testing <- l.pvalues
#     l.foldchanges_comparative.testing <- l.foldchanges_comparative
#     l.pvalues_comparative.testing <- l.pvalues_comparative
#   }
#   
#   # compute the actual values instead of fold changes
#   # individual comparative enrichment tests
#   
# }
# ### 
# 

print(lapply(l.foldchanges.testing[[m]][[1]], mean))
print(lapply(l.foldchanges.testing[[m]][[1]], mean))
print(lapply(l.foldchanges.training[[m]][[1]], mean))

for(m in 1:length(l.methods)){
  
  print(names(l.methods)[m])
  
  #l.foldchanges.training
  l.foldchanges_comparative.testing[[m]][[2]] <- lapply(l.foldchanges_comparative.testing[[m]][[2]], function(m) m[m != "Inf"])
  
  print(lapply(l.foldchanges_comparative.testing[[m]][[2]], mean))
  print(lapply(l.pvalues_comparative.testing[[m]][[2]], mean))
  
  # print(lapply(l.pvalues_comparative.testing[[m]][[2]], mean))
  
  #print(lapply(l.foldchanges.training[[m]][[2]], mean))
  
  # print(lapply(l.foldchanges.testing[[m]][[2]], mean))
  
  # lapply(l.foldchanges.training[[1]], mean)
  # lapply(l.foldchanges.training[[2]], mean)
}


# for supplement details and barplot
0.632 * unlist(lapply(l.foldchanges.testing[[1]], mean)) + 0.368 * unlist(lapply(l.foldchanges.training[[1]], mean)) # GRACE
0.632 * unlist(lapply(l.foldchanges.testing[[2]], mean)) + 0.368 * unlist(lapply(l.foldchanges.training[[2]], mean)) # initial network at similar size cutoff comparison (avoid naming GENIE3)

0.5 * unlist(lapply(l.pvalues.testing[[1]], mean)) + 0.5 * unlist(lapply(l.pvalues.training[[1]], mean))
0.5 * unlist(lapply(l.pvalues.testing[[2]], mean)) + 0.5 * unlist(lapply(l.pvalues.training[[2]], mean))

# GRACE vs ... 
# main barplots 





