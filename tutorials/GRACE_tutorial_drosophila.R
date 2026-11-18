
# GRACE 
# 
# Contact: Michael Banf, michael.banf@gmx.net	
# 
# This is the documentation for the GRACE algorithm. The implementation is a research prototype and is provided “as is”. No warranties or guarantees of any kind are provided. Do not distribute the GRACE R code or use it other than for your own research without permission by the author. 
# 
# GRACE has been written in R Version 3.2.2 - you might have to update to this version here https://cran.r-project.org


rm(list=ls()) # clear workspace 

# set path to the GRACE directory
setwd("/home/mbanf/Documents/Computational_Biology/Projects/GRACE_release/")

# general functions
source("GRACE.R")
source("GRACE_helperfunctions.R")

cat("\014")  
print("request and load libraries..")
request_libraries()

source("GRACE_optimization_routines.R")
source("GRACE_evaluation.R")
source("GRACE_tutorial_datasets.R")

### load tutorial datasets 
species = "drosophila"
print(paste("Species:", species))
lst.data <- request_dataset(species = species, method = "GENIE3") # a general tutorial

# regulators and targets
v.tfs.global <- lst.data$v.tfs
v.tgs.global <- lst.data$v.tgs

# initial high confidence network 
mat.gene_regulatory_network <- as.matrix(lst.data$mat.conserved_motif_plus_expression[v.tfs.global,v.tgs.global])
# mat.gene_regulatory_network <- mapping_intervals(mat.gene_regulatory_network, min(mat.gene_regulatory_network), max(mat.gene_regulatory_network), 0, 1) 

# co-functional network
mat.cofunctional_network <- lst.data$mat.FunctionalAssociation[v.tgs.global,v.tgs.global]

# Training, testing and validation sets`

# regulatory evidence based gold standard
mat.regulatory_evidence <- lst.data$mat.GS
#mat.regulatory_evidence <- lst.data$mat.GS
# co-functional evidence
mat.cofunctional_evidence <- lst.data$mat.GO.tgs[v.tgs.global,v.tgs.global]


lst.validation <- list(lst.data$mat.Chip, lst.data$mat.tissue_localization.tgs, lst.data$mat.HiC)
names(lst.validation) <- c("n.gs.grn.Chip","n.coreg.localization","n.coreg.HiC")
code.validation <- c("regulatory_evidence", "cofunctional_evidence", "cofunctional_evidence")


# additional validation sets (A. thaliana)
# mat.MR.subacon.tgs <- lst.data$mat.MR.subacon.tgs[v.tgs.global,v.tgs.global]
# mat.AtRegNet <- lst.data$mat.AtRegNet
# mat.aracyc.tgs <- lst.data$mat.aracyc.tgs[v.tgs.global,v.tgs.global]
# diag(mat.aracyc.tgs) <- 0

# global constraints
n.cpus <- 2
beta <- 1 # equal emphasis on precision (atrm) and recall (bp coreg pairs) - introducing (novel) combined f-measure (physical + functional evidence, singularity handling, minimum size)
th.min.coreg.tfs <- 1 # min number of TFs needed to define coregulation (jaccard > 0.5 leading to exclusive only for stringency) - obsolete
b.normalize_precision <- TRUE
b.jaccard = TRUE

n.sets <- 20
lambda.gridSearch <- c(0.01,seq(0.5,2.5,0.5))
#lambda.gridSearch <- seq(1.5,2,0.5)
#lambda.gridSearch <- seq(1.5,2.5,0.5)

# set.seed(1234)
# model learning related 
th.percentage.hyperparameters = 0.99 # grid search
max.call =  200 # simulated annealing 
n.models <- 100 # 100 models 
n.sample_size <- 0.632 # as in traditional bootrapping approaches

lst.benchmarks <- list(mat.regulatory_evidence, mat.cofunctional_evidence)
names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")

# precompute comparison model (=node potentials, solely based on gamma (v.th.grn))
if(FALSE){
  
  # needed for GRACE initiation - using random forest # 
  df.rate_density <- compute_fmeasures_regulatory_network(mat.grn=mat.gene_regulatory_network, 
                                                           lst.benchmarks=lst.benchmarks,
                                                           lst.validation=lst.validation,
                                                           code.validation = code.validation,
                                                           n.samples = 1000,
                                                           n.cpus=n.cpus) 
  
  # saveRDS(df.rate_density, "drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")

}else{
  df.rate_density <- readRDS("drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V3_unnormalized.rds")
}

#df.rate_density <- readRDS("drosophila_datasets/df.rate_density_GENIE3_jaccard_coreg_tfs_1_1000_V4_cell.0.rds")
lst.fscore_results <- plot_fmeasure_scores(df.rate_density)

# cairo_pdf("example.pdf", width = 7.5, height = 6, family= "Helvetica")#"DejaVu Sans")
lst.fscore_results$plot.PR 
# dev.off()

# constants needed globally
ratio.gs.grn_vs_bp.normalization <- lst.fscore_results$ratio.gs.grn_vs_bp.normalization
max.coreg.bp <- lst.fscore_results$max.coreg.bp
df.rate_density <- lst.fscore_results$df.rate_density

### variations ### 

# define limits 
gamma.upper <- alpha.upper <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.99999999)
#gamma.lower <- alpha.lower <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.9999)
####gamma.upper <- alpha.upper <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.9)
gamma.lower <- alpha.lower <- min(mat.gene_regulatory_network[mat.gene_regulatory_network>0]) #estimate_grn_cutoff(mat.grn=mat.grn.cns, v.cut=alpha.lower.cut)

gamma.gridSearch <- as.numeric(sapply(seq(1:n.sets), function(m) gamma.upper - (gamma.upper - gamma.lower)/n.sets * m))

th.preselection <- min(mat.gene_regulatory_network[mat.gene_regulatory_network>0]) #estimate_grn_cutoff(mat.grn=mat.grn.cns,v.cut=v.preselection)

# compute connectivity matrix per preslection (once)
print("construct meta co-regulation network")
lst.grace_support <- contruct_grace_meta_connectivity_network(mat.grn=mat.gene_regulatory_network, mat.FunctionalAssociation=mat.cofunctional_network, 
                                                              th.preselection = th.preselection, n.cores = n.cpus)

print("extract meta modules")
lst.modules <- prepare_meta_modules(adj.links=lst.grace_support$adj.links, df.grn = lst.grace_support$df.grn)

mat.grn.preselection <- acast(lst.grace_support$df.grn, TF~TG, value.var = "v.grn")
mat.grn.preselection[is.na(mat.grn.preselection)] <- 0
class(mat.grn.preselection) <- "numeric"   

v.tfs <- rownames(mat.grn.preselection)
v.tgs <- colnames(mat.grn.preselection)

mat.cofunctional_network <- mat.cofunctional_network[v.tgs,v.tgs]
mat.cofunctional_evidence <- mat.cofunctional_evidence[v.tgs,v.tgs]

####

tf <- intersect(rownames(mat.grn.preselection),rownames(mat.regulatory_evidence)) 
tgs <- intersect(colnames(mat.grn.preselection), colnames(mat.regulatory_evidence))

truereg.active <- mat.regulatory_evidence
truereg.active <- mat.grn.preselection[tf,tgs] * truereg.active[tf,tgs]
truereg.active[truereg.active < alpha.lower] <- 0
truereg.active[truereg.active > 0] <- 1
idx.tf <- names(which(rowSums(truereg.active) > 0))
idx.tg <- names(which(colSums(truereg.active) > 0))
truereg.active <- truereg.active[idx.tf,idx.tg]

truereg <- Matrix(0, nrow = length(rownames(mat.grn.preselection)), ncol = length(colnames(mat.grn.preselection)), 
                  dimnames = list(rownames(mat.grn.preselection), colnames(mat.grn.preselection)))
tfs <- intersect(rownames(truereg.active), rownames(truereg))
tgs <- intersect(colnames(truereg.active), colnames(truereg))
truereg[tfs, tgs] <- truereg.active[tfs, tgs]

idx.gs_grn <- which(truereg == 1)
idx.tg_bp <- which(mat.cofunctional_evidence == 1)

# 
# n.chip_background <- sum(lst.validation$n.gs.grn.Chip[v.tfs,v.tgs])
# n.chip_background <- sum(lst.validation$n.coreg.localization[v.tgs,v.tgs])
# n.chip_background <- sum(lst.validation$n.coreg.HiC[v.tgs,v.tgs])


if(FALSE){
#### multiple rounds of bootstrap aggregation - put into seperate method 
  strt<-Sys.time() 
  cl<-makeCluster(min(n.models, n.cpus))
  registerDoParallel(cl)
  
  ### running multiple bootstrap samples in parallel
  lst.grace_models <- foreach(i = 1:n.models, .packages=c("Matrix", "reshape2", "GenSA", "CRF", "igraph")) %dopar% {    
    
    source("GRACE_optimization_routines.R")
    source("GRACE.R")
    source("GRACE_helperfunctions.R")

    # use different seeds per run (simulated annealing)
    v.seed = 1234 + 158 * i
    
    lst.grace_models <- evaluate_mode_parameter(lst.grace_support=lst.grace_support, lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, 
                                                
                                                gamma.gridSearch=gamma.gridSearch,
                                                lambda.gridSearch=lambda.gridSearch,
                                                th.percentage.hyperparameters=th.percentage.hyperparameters,
                                                
                                                ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization,
                                                max.coreg.bp = max.coreg.bp,
                                                
                                                n.sample_size=n.sample_size, max.call=max.call, v.seed = v.seed, beta = beta)
    
    return(lst.grace_models)
    
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  #### remove models with fscore - 1 ###
  #saveRDS(lst.grace_models, "drosophila_datasets/lst.GRACE_models_genie3_m100_c200_beta1_V3.0_unnormalized.rds")
  #saveRDS(lst.grace_models, "drosophila_datasets/lst.GRACE_models_genie3_m100_c200_beta05_V3.0_unnormalized.rds")
  
}else{
  lst.grace_models <- readRDS("drosophila_datasets/lst.GRACE_models_genie3_m100_c200_beta1_V3.0_unnormalized.rds")
}  

lst.benchmarks <- list(mat.regulatory_evidence, mat.cofunctional_evidence)
names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")

###

print("Computing averaged bootstrap error validation - using same size network")
### compute the bootstrap averaged errors 


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
v.backgroundSet <- as.numeric(df.rate_density[nrow(df.rate_density),c(2,7,3,8,9)]) # needs to be adjusted
names(v.backgroundSet) <- validationSet

# generate training and test bootstrap prediction errors


n.links_background <- length(v.tfs) * length(v.tgs)
n.co_pairs <- length(v.tgs) * length(v.tgs)

n.gs_background <- sum(mat.regulatory_evidence[intersect(rownames(mat.regulatory_evidence), v.tfs), intersect(colnames(mat.regulatory_evidence), v.tgs)]) # 386 
n.chip_background <- sum(lst.validation$n.gs.grn.Chip[intersect(rownames(lst.validation$n.gs.grn.Chip), v.tfs), intersect(colnames(lst.validation$n.gs.grn.Chip), v.tgs)])

v.gs.grn <- c(n.gs_background, n.chip_background)
v.gs.copairs <- c(sum(mat.cofunctional_evidence[v.tgs,v.tgs]), sum(lst.validation$n.coreg.localization[v.tgs,v.tgs]), sum(lst.validation$n.coreg.HiC[v.tgs,v.tgs]))

####


l.foldchanges <- vector(mode = "list", length = 5)
l.pvalues <- vector(mode = "list", length = 5)

l.foldchanges[[1]] <- l.foldchanges[[2]] <- l.foldchanges[[3]] <- l.foldchanges[[4]] <- l.foldchanges[[5]] <- numeric(n.models)
l.pvalues[[1]] <- l.pvalues[[2]] <- l.pvalues[[3]] <- l.pvalues[[4]] <- l.pvalues[[5]] <- numeric(n.models)


b.validation <- TRUE


for(i in 1:n.models){
  
  cat("Processing... ", round(i/n.models * 100, digits = 2) , "%", "\r"); flush.console() 
  
  # use different seeds per run (simulated annealing)
  v.seed = 1234 + 158 * i
  set.seed(v.seed)
  
  #### Bootstrapping (based on regulators based meta module sampling (0.632))
  tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
  tgs.mrf <- unique(lst.grace_support$df.grn$TG)
  
  lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size, b.validation = b.validation)
  lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
  lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
  
  lst.benchmarks <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
  names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")

  lst.validation <- list(lst.data$mat.Chip, lst.data$mat.tissue_localization.tgs, lst.data$mat.HiC)
  names(lst.validation) <- c("n.gs.grn.Chip","n.coreg.localization","n.coreg.HiC")
  code.validation <- c("regulatory_evidence", "cofunctional_evidence", "cofunctional_evidence")
  
  ####
  
  alpha  <- lst.grace_models[[i]]$weights[1]
  gamma  <- lst.grace_models[[i]]$weights[2]
  lambda <- lst.grace_models[[i]]$weights[3]
  lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields(lst.grace_support=lst.grace_support, lst.modules=lst.modules, 
                                                            v.alpha = alpha, v.gamma=gamma, v.lambda=lambda, b.prim=TRUE) 
  mat.p.crf <- lst.mat.p.crf$mat.p.crf
  
  ### 
  
  
  l.eval <- compute_grace_complete_assessment(mat.grace=mat.p.crf, 
                                                lst.benchmarks=lst.benchmarks,
                                                lst.validation=lst.validation,
                                                code.validation = code.validation)


  n.links <- sum(mat.p.crf)
  

 # print(j)
  for(s in 1:length(v.sets_map)){
    
    hitInSample <- l.eval[[v.sets_map[s]]]
    
    if(s <= 2){
      hitInPop <- v.gs.grn[s]
      sampleSize <- n.links
      failInPop <- (n.links_background - hitInPop)
      foldchange <- (hitInSample/sampleSize)/(hitInPop/n.links_background)
    }else{
      hitInPop <- v.gs.copairs[s - 2]
      sampleSize <- l.eval$n.coreg.tg_pairs
      failInPop <- (n.co_pairs - hitInPop)
      foldchange <- (hitInSample/sampleSize)/(hitInPop/n.co_pairs)
    }
    
    #print(foldchange)
    
    p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    
    l.foldchanges[[s]][[i]] <- foldchange
    l.pvalues[[s]][[i]] <- p.val
    
  }
  #print("")
}


print("bootstrapping foldchange increases")
for(s in 1:length(v.sets_map)){
  print(v.sets_map[s])
  print(mean(l.foldchanges[[s]]))
}





print("build final ensemble model")
# model averaging - ensemble model building
alpha  <- lst.grace_models[[1]]$weights[1]
gamma  <- lst.grace_models[[1]]$weights[2]
lambda <- lst.grace_models[[1]]$weights[3]
lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields(lst.grace_support=lst.grace_support, lst.modules=lst.modules, 
                                                          v.alpha = alpha, v.gamma=gamma, v.lambda=lambda, b.prim=TRUE) 
mat.p.crf <- lst.mat.p.crf$mat.p.crf
n.models <- length(lst.grace_models)

for(j in 2:n.models){
  
  cat("Processing... ", round(j/n.models * 100, digits = 2) , "%", "\r"); flush.console()  
  alpha  <- lst.grace_models[[j]]$weights[1]
  gamma  <- lst.grace_models[[j]]$weights[2]
  lambda <- lst.grace_models[[j]]$weights[3]
  lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields(lst.grace_support=lst.grace_support, lst.modules=lst.modules, 
                                                            v.alpha = alpha, v.gamma=gamma, v.lambda=lambda, b.prim=TRUE) 
  mat.p.crf <- mat.p.crf + lst.mat.p.crf$mat.p.crf
}

mat.p.grace_ensemble <- mat.p.crf
mat.p.grace_ensemble <- mat.p.grace_ensemble/n.models

# define probabilistic cutoff {0.5 - 1}
# th.p.cutoff <- as.numeric(quantile(mat.p.grace_ensemble[mat.p.grace_ensemble>0], 0.5))
th.p.cutoff <- 0.5 # epsion parameter (defaul = 0.5, models = 100, simulated annealing runs)

if(th.p.cutoff == 0.5){
  mat.p.grace_ensemble[mat.p.grace_ensemble <= th.p.cutoff] <- 0
  mat.p.grace_ensemble[mat.p.grace_ensemble > th.p.cutoff] <- 1
}else{
  mat.p.grace_ensemble[mat.p.grace_ensemble < th.p.cutoff] <- 0
  mat.p.grace_ensemble[mat.p.grace_ensemble >= th.p.cutoff] <- 1
}



df.grace.ensemble <- as.data.frame(as.table(mat.p.grace_ensemble), stringsAsFactors = FALSE) 
names(df.grace.ensemble)[1:3] <- c("TF","TG","value")
df.grace.ensemble <- subset(df.grace.ensemble, df.grace.ensemble$value > 0)
df.grace.ensemble <- df.grace.ensemble[order(-df.grace.ensemble$value),] 

write.csv(df.grace.ensemble[,1:2], "prediction_drosophila_graceEnsemble.csv", row.names = FALSE)


