
# A. thaliana based tutorial of the GRACE application #
# GRACE version 1.0 (AraGRACE)
# 
# to be used with other Arabidopsis specific inferred gene regulatory network (for general datasets use GRACE version 1.1. (the drosophila tutorial). GRACE version 1.1 has been made faster.. )
# 
# a) call infer network to predict gene regulatory network (parameters: gene expression + binding dataframe)
# b) then replace "mat.gene_regulatory_network" with inferred network and run the rest of the pipeline

rm(list=ls()) # clear workspace 

# set path to the GRACE directory
setwd("/home/mbanf/Documents/Computational_Biology/Projects/GRACE_release/")

source("GRACE.R")
source("GRACE_helperfunctions.R")

cat("\014")  
print("request and load libraries..")
request_libraries()

source("GRACE_evaluation.R")
source("GRACE_tutorial_datasets.R")


### load tutorial datasets 
species = "arabidopsis"
print(paste("Species:", species))
lst.data <- request_dataset(species = species)

# regulators and targets
v.tfs.global <- lst.data$v.tfs
v.tgs.global <- lst.data$v.tgs

# initial high confidence network 
mat.gene_regulatory_network <- lst.data$mat.conserved_motif_plus_expression[v.tfs.global,v.tgs.global]

# co-functional network
mat.cofunctional_network <- lst.data$mat.FunctionalAssociation[v.tgs.global,v.tgs.global]

# Training, testing and validation sets`

# regulatory evidence based gold standard
mat.regulatory_evidence <- lst.data$mat.GS
# co-functional evidence
mat.cofunctional_evidence <- lst.data$mat.GO.tgs[v.tgs.global,v.tgs.global]

# additional validation sets (A. thaliana)
mat.MR.subacon.tgs <- lst.data$mat.MR.subacon.tgs[v.tgs.global,v.tgs.global]
mat.AtRegNet <- lst.data$mat.AtRegNet
mat.aracyc.tgs <- lst.data$mat.aracyc.tgs[v.tgs.global,v.tgs.global]
diag(mat.aracyc.tgs) <- 0


# global constraints
n.cpus <- 2 #14
beta <- 1 # equal emphasis on precision (atrm) and recall (bp coreg pairs) - introducing (novel) combined f-measure (physical + functional evidence, singularity handling, minimum size)
n.models <- 100 # 100 models 
n.sample_size <- 0.632 # as in traditional bootrapping approaches

th.min.coreg.tfs <- 1 # min number of TFs needed to define coregulation (jaccard > 0.5 leading to exclusive only for stringency)
b.normalize_precision <- TRUE
b.jaccard = TRUE

n.sets <- 20
lambda.gridSearch <- c(0.01,seq(0.5,2.5,0.5))


set.seed(1234)
# model learning related 
th.percentage.hyperparameters = 0.99 # grid search
max.call =  250 # simulated annealing

lst.benchmarks <- list(mat.regulatory_evidence, mat.cofunctional_evidence)
names(lst.benchmarks) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")

# precomputed comparison model (=node potentials, solely based on gamma (v.th.grn))
if(FALSE){
  df.rate_density <- compute_fmeasures_arabidopsis(mat.grn=mat.gene_regulatory_network, 
                                                   lst.benchmarks=lst.benchmarks,
                                                   
                                                   mat.AtRegNet=mat.AtRegNet,
                                                   mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                   mat.aracyc.tgs=mat.aracyc.tgs,
                                                   
                                                   b.jaccard = b.jaccard,
                                                   th.min.coreg.tfs = th.min.coreg.tfs, 
                                                   n.samples = 100,
                                                   n.cpus=n.cpus) 
  
  #saveRDS(df.rate_density, "athaliana_datasets/df.rate_density_de_jaccard_coreg_tfs_1_1000_032316.rds")
  
}else{
  df.rate_density <- readRDS("athaliana_datasets/df.rate_density_de_jaccard_coreg_tfs_1_1000.rds")

}

max.gs.grn <- max(df.rate_density$n.gs.grn)
max.coreg.tgs.bp <- max(df.rate_density$n.coreg.tg_pairs.BP)

max(df.rate_density$n.grn_agris.gs)
max(df.rate_density$n.tg_pairs.bp)
max(df.rate_density$n.coreg.tg_pairs.suba)
max(df.rate_density$n.coreg.tg_pairs.aracyc)

# extract ratio of precision x recall
Precision.gs.grn <- (df.rate_density$n.gs.grn/df.rate_density$n.pred)
#Recall.coreg.tgs.bp <- df.rate_density$n.coreg.tg_pairs.BP/max(df.rate_density$n.coreg.tg_pairs.BP)

ratio.gs.grn_vs_bp.normalization <- 1 / max(Precision.gs.grn) #/ max(Recall.coreg.tgs.bp)
max.coreg.bp <- max(df.rate_density$n.coreg.tg_pairs.BP)

####
x <- df.rate_density$n.pred / max(df.rate_density$n.pred)

y1 <- (df.rate_density$n.gs.grn/df.rate_density$n.pred)
y1[is.na(y1)] <- 0
y1 <- y1/max(y1) # normalization
#plot(x,y1, type = "l") # training (ATRM)

y2 <- df.rate_density$n.coreg.tg_pairs.BP/max(df.rate_density$n.coreg.tg_pairs.BP)
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



# cairo_pdf("example.pdf", width = 7.5, height = 6, family= "Helvetica")#"DejaVu Sans")
ggplot(data, aes(x=x, y=y)) +
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
            theme(legend.justification=c(1,0), legend.position=c(1,0))
# dev.off()


# define limits 
gamma.upper <- alpha.upper <- quantile(mat.gene_regulatory_network[mat.gene_regulatory_network>0], 0.99999999)
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

mat.cofunctional_network <- lst.data$mat.FunctionalAssociation[v.tgs,v.tgs]
mat.cofunctional_evidence <- lst.data$mat.GO.tgs[v.tgs,v.tgs]

# validation sets
mat.MR.subacon.tgs <- lst.data$mat.MR.subacon.tgs[v.tgs,v.tgs]
mat.AtRegNet <- lst.data$mat.AtRegNet
mat.aracyc.tgs <- lst.data$mat.aracyc.tgs[v.tgs,v.tgs]
diag(mat.aracyc.tgs) <- 0

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


if(FALSE){ 
#### multiple rounds of bootstrap aggregation - put into seperate method 
  strt<-Sys.time() 
  cl<-makeCluster(min(n.models, n.cpus))
  registerDoParallel(cl)
  
  ### running multiple bootstrap samples in parallel
  lst.grace_models <- foreach(i = 1:n.models, .packages=c("Matrix", "reshape2", "GenSA", "CRF", "igraph")) %dopar% {    
    
    #for(i in 1:n.models){
    
    # optimiation function - grid search
    # optimiation function - grid search
    gridSearchOptimization <- function(gamma.gridSearch, lambda.gridSearch, th.percentage.hyperparameters = 0.99){
      
      df.gridSearch <- expand.grid(gamma.gridSearch, lambda.gridSearch)
      names(df.gridSearch) <- c("gamma", "lambda")
      
      df.gridSearch["f_score"] <- 0
      df.gridSearch["n.predictions"] <- 0
      df.gridSearch["n.gs.grn"] <- 0
      df.gridSearch["n.coreg.bp.tgs"] <- 0
      
      for(j in 1:nrow(df.gridSearch)){
        
        cat("Processing... ", round(j/nrow(df.gridSearch) * 100, digits = 2) , "%", "\r"); flush.console()  
        v.alpha = v.gamma = df.gridSearch$gamma[j]; v.lambda= df.gridSearch$lambda[j];
        
        if(v.alpha >= v.gamma){
          lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                                     v.alpha = v.alpha, v.gamma = v.gamma, v.lambda = v.lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
          
          mat.p.crf <- lst.mat.p.crf$mat.p.crf
          lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
                                                      lst.benchmarks=lst.gs.grn_bp.train,
                                                      th.min.coreg.tfs = th.min.coreg.tfs, 
                                                      b.jaccard = b.jaccard,
                                                      beta = beta)
          
          y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
          y1[is.na(y1)] <- 0
          if(b.normalize_precision){
            y1 <- ratio.gs.grn_vs_bp.normalization * y1  
          }
          
          y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
          y2[is.na(y2)] <- 0
          
          fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
          fscore_beta[is.na(fscore_beta)] <- 0
          
          fscore <- fscore_beta
          n.gs.grn <- lst.eval$n.gs.grn
          n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
          n.pred <- lst.eval$n.pred
          
        }else{
          fscore <- - 100
          n.gs.grn <- - 1
          n.coreg.bp.tgs <- - 1
          n.pred <- - 1
        }
        
        df.gridSearch$f_score[j] <- fscore
        df.gridSearch$n.predictions[j] <- n.pred
        df.gridSearch$n.gs.grn[j] <- n.gs.grn
        df.gridSearch$n.coreg.bp.tgs[j] <- n.coreg.bp.tgs
        #print(paste("alpha", v.alpha, ", gamma", v.gamma, ", lambda", v.lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
      }
      
      # return Simulated annealing range 
      v.fscore.selection <- quantile(df.gridSearch$f_score, th.percentage.hyperparameters)
      df.gridSearch <- subset(df.gridSearch, df.gridSearch$f_score >= v.fscore.selection)
      
      v.double_fscores <- names(which(table(df.gridSearch$f_score) > 1))
      if(length(v.double_fscores) > 0){
        for(l in 1:length(v.double_fscores)){
          v.lambda.tmp <- subset(df.gridSearch, df.gridSearch$f_score == v.double_fscores[l])$lambda
          v.lambda.tmp <- v.lambda.tmp[!v.lambda.tmp %in% min(v.lambda.tmp)]
          
          idx.remove <- which(df.gridSearch$f_score == v.double_fscores[l] & df.gridSearch$lambda %in% v.lambda.tmp)
          df.gridSearch <- df.gridSearch[-idx.remove,]
        }
      }
      
      v.gamma.selection <- unique(df.gridSearch$gamma)
      v.lambda.selection <- unique(df.gridSearch$lambda)
      
      w  <- numeric(2)
      lower_limits <- numeric(2)
      upper_limits <- numeric(2)
      
      # gamma - single value handling
      if(length(v.gamma.selection) == 1){
        idx.gamma <- which(gamma.gridSearch == v.gamma.selection)
        if(idx.gamma == 1){
          lower_limits[[1]] <- v.gamma.selection
          upper_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma - 1]) # inverse indexing because of v.grn ordering
        }else if(idx.gamma == length(gamma.gridSearch)){
          lower_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma + 1])
          upper_limits[[1]] <- v.gamma.selection
        }else{
          lower_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma + 1])
          upper_limits[[1]] <- as.numeric(gamma.gridSearch[idx.gamma - 1])
        }
        w[[1]] <- v.gamma.selection
      }else{
        w[[1]] <- min(v.gamma.selection)
        lower_limits[[1]] <- min(v.gamma.selection)
        upper_limits[[1]] <- max(v.gamma.selection)
      }
      
      # lambda - single value handling
      if(length(v.lambda.selection) == 1){
        idx.lambda <- which(lambda.gridSearch == v.lambda.selection)
        if(idx.lambda == 1){
          lower_limits[[2]] <- v.lambda.selection
          upper_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda + 1])
        }else if(idx.lambda == length(lambda.gridSearch)){
          lower_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda - 1])
          upper_limits[[2]] <- v.lambda.selection
        }else{
          lower_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda - 1])
          upper_limits[[2]] <- as.numeric(lambda.gridSearch[idx.lambda + 1])
        }
        w[[2]] <- v.lambda.selection
      }else{
        w[[2]] <- min(v.lambda.selection)
        lower_limits[[2]] <- min(v.lambda.selection)
        upper_limits[[2]] <- max(v.lambda.selection)
      }
      
      return(list(w=w,lower_limits=lower_limits,upper_limits=upper_limits))
    }
    
    # optimization function - Simulated annealing
    evaluate_SA <- function(w) {
      
      v.alpha = v.gamma = w[1]; v.lambda= w[2]; 
      if(v.alpha >= v.gamma){
        lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                                   v.alpha = v.alpha, v.gamma = v.gamma, v.lambda = v.lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
        
        mat.p.crf <- lst.mat.p.crf$mat.p.crf
        lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
                                                    lst.benchmarks=lst.gs.grn_bp.train,
                                                    th.min.coreg.tfs = th.min.coreg.tfs, 
                                                    b.jaccard = b.jaccard,
                                                    beta = beta)
        
        y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
        y1[is.na(y1)] <- 0
        if(b.normalize_precision){
          y1 <- ratio.gs.grn_vs_bp.normalization * y1  
        }
        
        y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
        y2[is.na(y2)] <- 0
        
        fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
        fscore_beta[is.na(fscore_beta)] <- 0
        
        fscore <- fscore_beta
        n.gs.grn <- lst.eval$n.gs.grn
        n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
        n.pred <- lst.eval$n.pred
        
      }else{
        fscore <- - 100
        n.gs.grn <- - 1
        n.coreg.bp.tgs <- - 1
        n.pred <- - 1
      }
      
      print(paste("alpha", v.alpha, ", gamma", v.gamma, ", lambda", v.lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
      returnVal = fscore * -1 
      returnVal
    }
    
    # use different seeds per run (simulated annealing)
    v.seed = 1234 + 158 * i
    set.seed(v.seed)
    
    #### Bootstrapping (based on regulators based meta module sampling (0.632))
    
    tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
    tgs.mrf <- unique(lst.grace_support$df.grn$TG)
    
    lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size)
    lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
    lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
    
    lst.gs.grn_bp.train <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
    names(lst.gs.grn_bp.train) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    ## Step 1: Grid Search (Coarse) optimization to identify finer range for Simulated annealing optimization ##
    strt<-Sys.time() 
    lst.hyperparameter <- gridSearchOptimization(gamma.gridSearch, lambda.gridSearch, th.percentage.hyperparameters)
    w.intials     <- lst.hyperparameter$w
    lower_limits <- lst.hyperparameter$lower_limits
    upper_limits <- lst.hyperparameter$upper_limits
    print(Sys.time()-strt)
    #     w     <- c(gamma.upper, 1.5) #0.5) # start at upper end - therefore, the one with smallest 
    #     lower <- c(gamma.lower, 2.0) #0.01) 
    #     upper <- c(gamma.upper, 2.5) # 1.0)
    
    ## Step2: Simulated annealing ##
    strt<-Sys.time() 
  
    out <- GenSA(w.intials, lower = lower_limits, upper = upper_limits, fn = evaluate_SA , control=list(max.call=max.call)) 
    #out <- GenSA(w, lower = lower, upper = upper, fn = evaluate_SA , control=list(max.call=max.call, smooth = FALSE, simple.function = TRUE))
    results <- out[c("value", "par", "counts")]
    weights <- results$par # go for smallest lambda
    
    print(Sys.time()-strt)
    
    #### construct model using learned weights
    alpha <- gamma <- weights[1]
    lambda <- weights[2]
    
    lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                               v.alpha = alpha, v.gamma = gamma, v.lambda = lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
    
    mat.p.crf <- lst.mat.p.crf$mat.p.crf
    lst.gs.grn_bp <- list(truereg, mat.cofunctional_evidence)
    names(lst.gs.grn_bp) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
                                                lst.benchmarks=lst.gs.grn_bp,  
                                                th.min.coreg.tfs = th.min.coreg.tfs, 
                                                b.jaccard = b.jaccard,
                                                beta = beta)
    
    y1 <- (lst.eval$n.gs.grn/lst.eval$n.pred)
    y1[is.na(y1)] <- 0
    if(b.normalize_precision){
      y1 <- ratio.gs.grn_vs_bp.normalization * y1  
    }
    
    y2 <- lst.eval$n.coreg.bp_tg_pairs/max.coreg.bp
    y2[is.na(y2)] <- 0
    
    fscore_beta <- (1 + beta^2)*(y1*y2/((beta^2 * y1)+y2))
    fscore_beta[is.na(fscore_beta)] <- 0
    
    fscore <- fscore_beta
    n.gs.grn <- lst.eval$n.gs.grn
    n.coreg.bp.tgs <- lst.eval$n.coreg.bp_tg_pairs
    n.pred <- lst.eval$n.pred
    
    print(paste("alpha", alpha, ", gamma", gamma, ", lambda", lambda, ", f_score", fscore, "n.predictions", n.pred, "n.gs.grn",n.gs.grn ,"n.coreg.bp.tgs", n.coreg.bp.tgs))
    
    lst.grace_models <- list(weights=c(alpha, gamma, lambda), fscore=fscore, n.pred=n.pred,n.gs.grn=n.gs.grn,n.coreg.bp.tgs=n.coreg.bp.tgs)
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
}else{
  # individual learned models used to generate figure 4 in publication 
  lst.grace_models <- readRDS("athaliana_datasets/lst.GRACE_models_100_100_250_032316.rds")
}  


print("Computing averaged bootstrap error validation - using same size network")
### compute the bootstrap averaged errors 

l.foldchanges <- vector(mode = "list", length = 2)
names(l.foldchanges) <- c("GRACE", "Original")

for(j in 1:2){
  l.foldchanges[[j]] <- vector(mode = "list", length = 5)
  for(s in 1:5){
    l.foldchanges[[j]][[s]] <- numeric(n.models)
  }
}

l.pvalues <- vector(mode = "list", length = 2)
names(l.pvalues) <- c("GRACE", "Original")

for(j in 1:2){
  l.pvalues[[j]] <- vector(mode = "list", length = 5)
  for(s in 1:5){
    l.pvalues[[j]][[s]] <- numeric(n.models)
  }
}


validationSet = c("Gene regulatory links\n (AGRIS)", 
                  "Gene regulatory links\n (ATRM)", 
                  "Co-functional pairs\n (Gene Ontology)",
                  "Co-localized pairs\n (Suba)", 
                  "Metabolic pathway pairs\n (Aracyc)")

v.sets_map <- c(1,3,4,8,9)
names(v.sets_map) <- validationSet

v.background = df.rate_density[nrow(df.rate_density),]
n.pred.background <- df.rate_density$n.pred[nrow(df.rate_density)]
n.coreg.background <- df.rate_density$n.coreg.tg_pairs[nrow(df.rate_density)]
v.backgroundSet <- as.numeric(df.rate_density[nrow(df.rate_density),c(1,6,7,15,16)]) # needs to be adjusted
names(v.backgroundSet) <- validationSet


#n.links_background <- length(v.tfs) * length(v.tgs)
#n.co_pairs <- length(v.tgs) * length(v.tgs)

n.gs_background <- sum(mat.regulatory_evidence[intersect(rownames(mat.regulatory_evidence), v.tfs), intersect(colnames(mat.regulatory_evidence), v.tgs)]) # 386 
n.chip_background <- sum(mat.AtRegNet[intersect(rownames(mat.AtRegNet), v.tfs), intersect(colnames(mat.AtRegNet), v.tgs)])

v.gs.grn <- c(n.chip_background, n.gs_background)
v.gs.copairs <- c(sum(lst.data$mat.GO.tgs[v.tgs,v.tgs]), sum(lst.data$mat.MR.subacon.tgs[v.tgs,v.tgs]), sum(lst.data$mat.aracyc.tgs[v.tgs,v.tgs]))


n.links_background <- length(v.tfs) * length(v.tgs)
n.co_pairs <- length(v.tgs) * length(v.tgs)



# generate training and test bootstrap prediction errors
for(t in 1:2){
  
  if(t == 1){
    b.validation <- FALSE
  }else{
    b.validation <- TRUE
  }

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
    
    lst.validation <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
    names(lst.validation) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    ####
    
    alpha  <- lst.grace_models[[i]]$weights[1]
    gamma  <- lst.grace_models[[i]]$weights[2]
    lambda <- lst.grace_models[[i]]$weights[3]
    lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields(lst.grace_support=lst.grace_support, lst.modules=lst.modules, 
                                                              v.alpha = alpha, v.gamma=gamma, v.lambda=lambda, b.prim=TRUE) 
    mat.p.crf <- lst.mat.p.crf$mat.p.crf
    
    
    ###
     
    l.eval <- vector(mode = "list", length = 2)
    
    l.eval[[1]] <-  compute_grace_ensemble_model_assessment(mat.p.grace_ensemble=mat.p.crf, 
                                                         
                                                             lst.benchmarks=lst.validation, # complementary validation set per bootstrap model
                                                             
                                                             mat.AtRegNet=mat.AtRegNet,
                                                             mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                             mat.aracyc.tgs=mat.aracyc.tgs,
                                                             
                                                             b.jaccard = TRUE,
                                                             th.min.coreg.tfs = 1, 
                                                             beta = beta)
        
        
    n.links <- sum(mat.p.crf)
    
    
    
    
    
    for(m in 1:length(l.methods)){
      
      mat.gene_regulatory_network <- as.matrix(l.methods[[m]])
    
      mat.grn.comparison <- mat.gene_regulatory_network[v.tfs,v.tgs]
      # mat.grn.comparison <- mat.grn.comparison[rownames(lst.validation$mat.regulatory_evidence),]
      
      df.grn.comparison <- as.data.frame(as.table(mat.grn.comparison), stringsAsFactors = FALSE) 
      names(df.grn.comparison)[1:3] <- c("TF","TG","v.grn")
      df.grn.comparison <- subset(df.grn.comparison, df.grn.comparison$v.grn > 0)
      df.grn.comparison <- df.grn.comparison[order(-df.grn.comparison$v.grn),] 
      df.grn.comparison <- df.grn.comparison[1:n.links,]
      
      mat.grn.comparison <- acast(df.grn.comparison, TF~TG, value.var="v.grn")
      mat.grn.comparison[is.na(mat.grn.comparison)] <- 0; 
      mat.grn.comparison[mat.grn.comparison > 0] <- 1
      
      
      l.eval[[2]] <-  compute_grace_ensemble_model_assessment(mat.p.grace_ensemble=mat.grn.comparison, 
                                                           
                                                           lst.benchmarks=lst.validation, # complementary test set set per bootstrap model
                                                           
                                                           mat.AtRegNet=mat.AtRegNet, # complete validation sets
                                                           mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                           mat.aracyc.tgs=mat.aracyc.tgs,
                                                           
                                                           b.jaccard = TRUE,
                                                           th.min.coreg.tfs = 1, 
                                                           beta = beta)
      
      for(j in 1:2){
        
        for(s in 1:length(v.sets_map)){
          
          hitInSample <- l.eval[[j]][[v.sets_map[s]]]
          
          if(s <= 2){
            hitInPop <- v.gs.grn[s]
            sampleSize <- n.links
            failInPop <- (n.links_background - hitInPop)
            foldchange <- (hitInSample/sampleSize)/(hitInPop/n.links_background)
          }else{
            hitInPop <- v.gs.copairs[s - 2]
            sampleSize <- l.eval[[j]]$n.coreg.tg_pairs
            failInPop <- (n.co_pairs - hitInPop)
            foldchange <- (hitInSample/sampleSize)/(hitInPop/n.co_pairs)
          }
          
          #print(foldchange)
          
          p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
          
          l.foldchanges[[m]][[j]][[s]][[i]] <- foldchange
          l.pvalues[[m]][[j]][[s]][[i]] <- p.val
          
        }
      }
      
      for(s in 1:length(v.sets_map)){
        
        hitInSample <- l.eval[[1]][[v.sets_map[s]]]
        hitInPop <- l.eval[[2]][[v.sets_map[s]]]
        
        if(s <= 2){
          sampleSize <- n.links
          failInPop <- (n.links - hitInPop)
          foldchange <- (hitInSample/n.links)/(hitInPop/n.links)
        }else{
          sampleSize <- l.eval[[1]]$n.coreg.tg_pairs
          failInPop <- (l.eval[[2]]$n.coreg.tg_pairs - hitInPop)
          foldchange <- (hitInSample/sampleSize)/(hitInPop/l.eval[[2]]$n.coreg.tg_pairs)
        }
        
        tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
        p.val <- fisher.test(tab)$p.value
        
        l.foldchanges_comparative[[m]][[j]][[s]][[i]] <- foldchange
        l.pvalues_comparative[[m]][[j]][[s]][[i]] <- p.val
        
      }
      
    }
  }
    
  if(t == 1){
    l.foldchanges.training <- l.foldchanges
    l.pvalues.training <- l.pvalues
    l.foldchanges_comparative.training <- l.foldchanges_comparative
    l.pvalues_comparative.training <- l.pvalues_comparative
  }else{
    l.foldchanges.testing <- l.foldchanges
    l.pvalues.testing <- l.pvalues
    l.foldchanges_comparative.testing <- l.foldchanges_comparative
    l.pvalues_comparative.testing <- l.pvalues_comparative
  }

}
### 



lapply(l.pvalues.training[[2]], mean)
lapply(l.pvalues.testing[[2]], mean)

lapply(l.foldchanges.training[[1]], mean)
lapply(l.foldchanges.training[[2]], mean)


# for supplement details and barplot
0.632 * unlist(lapply(l.foldchanges.testing[[1]], mean)) + 0.368 * unlist(lapply(l.foldchanges.training[[1]], mean)) # GRACE
0.632 * unlist(lapply(l.foldchanges.testing[[2]], mean)) + 0.368 * unlist(lapply(l.foldchanges.training[[2]], mean)) # initial network at similar size cutoff comparison (avoid naming GENIE3)

0.5 * unlist(lapply(l.pvalues.testing[[1]], mean)) + 0.5 * unlist(lapply(l.pvalues.training[[1]], mean))
0.5 * unlist(lapply(l.pvalues.testing[[2]], mean)) + 0.5 * unlist(lapply(l.pvalues.training[[2]], mean))

# GRACE vs ... 
# main barplots 










lapply(l.foldchanges.testing[[1]], function(m) mean(unlist(m)))
lapply(l.foldchanges.testing[[2]], function(m) mean(unlist(m)))




# significant increase?


# df.validation["asterisk"] = ifelse(df.validation$p.val < 0.001, "**", ifelse(df.validation$p.val < 0.01, "*", ""))#, ifelse(df.validation$p.val < 0.05, "*", "")))


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

# relative improvemnet 


# 
#     
#     
#     for(k in 1:length(v.validation.node_potentals)){
#       
#       hitInSample <- lst.res.groups[[k]][j]
#       sampleSize <- lst.res.groups.backgrounds[[k]][j]
#       
#       hitInPop <- v.validation.node_potentals[j]
#       failInPop <- (v.validation_background.node_potentials[j] - hitInPop)
#       
#       foldchange <- (hitInSample/sampleSize)/(hitInPop/v.validation_background.node_potentials[j])
#       
#       if(m.test == "hyper"){
#         p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
#       }
#       
#       df.validation <- rbind(df.validation, data.frame(group = v.groups[k],
#                                                        validationSet = validationSet[j],
#                                                        foldChange = foldchange,
#                                                        p.value = p.val))
#       
#     }
#     
# 
#     
#     v.groups = c(paste("GRACE's link prioritizations (",v.network_sizes[3]," links)", sep = ""),
#                  paste("High confidence network 1 (GENIE3 + conserved motifs, ",v.network_sizes[2]," links)", sep = ""),
#                  paste("High confidence network 2 (GENIE3 + conserved motifs, ",v.network_sizes[1]," links)", sep = ""))
#     
#     #   v.groups = c(paste("GRACE's link prioritizations ( ",v.network_sizes[3]," links)", sep = ""),
#     #                paste("High confidence network 1 ( ",v.network_sizes[2]," links)", sep = ""),
#     #                paste("High confidence network 2 ( ",v.network_sizes[1]," links)", sep = ""))
#     
#     lst.res.groups <- vector(mode = "list", length = 3)
#     lst.res.groups.backgrounds <- vector(mode = "list", length = 3)
#     
#     lst.res.groups[[1]] <- v.validation.grace
#     lst.res.groups[[2]] <- v.validation.node_potentials.equal_size
#     lst.res.groups[[3]] <- v.validation.node_potentials.best_fscore
#     
#     lst.res.groups.backgrounds[[1]] <- v.validation_background.grace
#     lst.res.groups.backgrounds[[2]] <- v.validation_background.node_potentials.equal_size
#     lst.res.groups.backgrounds[[3]] <- v.validation_background.node_potentials.best_fscore
#     
#     
#     #   v.groups = c(paste("GRACE (", round(lst.eval_grace$n.pred), " links)", sep = ""),
#     #                paste("conserved motifs + GENIE3 (", round(lst.eval_grace$n.pred), " links)", sep = ""),
#     #                paste("conserved motifs + GENIE3 (", round(lst.eval_grace$n.pred), " links)", sep = ""))
#     
#     
#     
#     df.validation <- data.frame(group = character(),
#                                 validationSet = character(),
#                                 foldChange = numeric(),
#                                 p.value = numeric())
#     
#     for(k in length(v.groups):1){
#       for(j in 1:length(v.validation.node_potentals)){
#         
#         hitInSample <- lst.res.groups[[k]][j]
#         sampleSize <- lst.res.groups.backgrounds[[k]][j]
#         
#         hitInPop <- v.validation.node_potentals[j]
#         failInPop <- (v.validation_background.node_potentials[j] - hitInPop)
#         
#         foldchange <- (hitInSample/sampleSize)/(hitInPop/v.validation_background.node_potentials[j])
#         
#         if(m.test == "hyper"){
#           p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
#         }
#         
#         df.validation <- rbind(df.validation, data.frame(group = v.groups[k],
#                                                          validationSet = validationSet[j],
#                                                          foldChange = foldchange,
#                                                          p.value = p.val))
#         
#       }
#     }
#     
#     df.validation["asterisk"] = ifelse(df.validation$p.val < 0.001, "**", ifelse(df.validation$p.val < 0.01, "*", ""))#, ifelse(df.validation$p.val < 0.05, "*", "")))
#     
#     
#     
#     
#     
#     
#     vector for the results to be averaged (with error bars... )
#     
#     # fold change enrichment over all methods
#     v.FCs
#     
#     
#   }
#   
# }

# average fold change enrichment of improvments






print("Constructing final ensemble model")

#   lst.gs.grn_bp <- list(truereg, mat.tgs.aranet_benchmark)
#   names(lst.gs.grn_bp) <- c("mat.GS", "mat.BP.coreg")

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

write.csv(df.grace.ensemble[,1:2], "prediction_arabidopsis_graceEnsemble.csv", row.names = FALSE)



# evaluate the ensemble model
l.eval <- vector(mode = "list", length = 2)

l.eval[[1]]  <-  compute_grace_ensemble_model_assessment(mat.p.grace_ensemble=mat.p.grace_ensemble, 
                                                         
                                                         lst.benchmarks=lst.benchmarks,
                                                         
                                                         mat.AtRegNet=mat.AtRegNet,
                                                         mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                         mat.aracyc.tgs=mat.aracyc.tgs,
                                                         
                                                         b.jaccard = TRUE,
                                                         th.min.coreg.tfs = 1, 
                                                         beta = beta)

  
n.links <- sum(mat.p.grace_ensemble)

for(m in 1:length(l.methods)){
  
  mat.gene_regulatory_network <- as.matrix(l.methods[[m]])
  
  mat.grn.comparison <- mat.gene_regulatory_network[v.tfs,v.tgs]
  # mat.grn.comparison <- mat.grn.comparison[rownames(lst.validation$mat.regulatory_evidence),]
  
  df.grn.comparison <- as.data.frame(as.table(mat.grn.comparison), stringsAsFactors = FALSE) 
  names(df.grn.comparison)[1:3] <- c("TF","TG","v.grn")
  df.grn.comparison <- subset(df.grn.comparison, df.grn.comparison$v.grn > 0)
  df.grn.comparison <- df.grn.comparison[order(-df.grn.comparison$v.grn),] 
  df.grn.comparison <- df.grn.comparison[1:n.links,]
  
  mat.grn.comparison <- acast(df.grn.comparison, TF~TG, value.var="v.grn")
  mat.grn.comparison[is.na(mat.grn.comparison)] <- 0; 
  mat.grn.comparison[mat.grn.comparison > 0] <- 1
  
  l.eval[[2]] <-  compute_grace_ensemble_model_assessment(mat.p.grace_ensemble=mat.grn.comparison, 
                                                          
                                                          lst.benchmarks=lst.benchmarks, # complementary test set set per bootstrap model
                                                          
                                                          mat.AtRegNet=mat.AtRegNet, # complete validation sets
                                                          mat.MR.subacon.tgs=mat.MR.subacon.tgs,
                                                          mat.aracyc.tgs=mat.aracyc.tgs,
                                                          
                                                          b.jaccard = TRUE,
                                                          th.min.coreg.tfs = 1, 
                                                          beta = beta)
  
  for(j in 1:2){
    
    for(s in 1:length(v.sets_map)){
      
      hitInSample <- l.eval[[j]][[v.sets_map[s]]]
      
      if(s <= 2){
        hitInPop <- v.gs.grn[s]
        sampleSize <- n.links
        failInPop <- (n.links_background - hitInPop)
        foldchange <- (hitInSample/sampleSize)/(hitInPop/n.links_background)
      }else{
        hitInPop <- v.gs.copairs[s - 2]
        sampleSize <- l.eval[[j]]$n.coreg.tg_pairs
        failInPop <- (n.co_pairs - hitInPop)
        foldchange <- (hitInSample/sampleSize)/(hitInPop/n.co_pairs)
      }
      
      #print(foldchange)
      
      p.val <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
      
      l.foldchanges[[m]][[j]][[s]][[i]] <- foldchange
      l.pvalues[[m]][[j]][[s]][[i]] <- p.val
      
    }
  }
  
  for(s in 1:length(v.sets_map)){
    
    hitInSample <- l.eval[[1]][[v.sets_map[s]]]
    hitInPop <- l.eval[[2]][[v.sets_map[s]]]
    
    if(s <= 2){
      sampleSize <- n.links
      failInPop <- (n.links - hitInPop)
      foldchange <- (hitInSample/n.links)/(hitInPop/n.links)
    }else{
      sampleSize <- l.eval[[1]]$n.coreg.tg_pairs
      failInPop <- (l.eval[[2]]$n.coreg.tg_pairs - hitInPop)
      foldchange <- (hitInSample/sampleSize)/(hitInPop/l.eval[[2]]$n.coreg.tg_pairs)
    }
    
    tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, failInPop), nrow = 2)
    p.val <- fisher.test(tab)$p.value
    
    l.foldchanges_comparative[[m]][[j]][[s]][[i]] <- foldchange
    l.pvalues_comparative[[m]][[j]][[s]][[i]] <- p.val
    
  }
  
}


l.foldchanges.testing <- l.foldchanges
l.pvalues.testing <- l.pvalues
l.foldchanges_comparative.testing <- l.foldchanges_comparative
l.pvalues_comparative.testing <- l.pvalues_comparative
### 






# 
# # that is the comparison model - equal size, f score are not really comparable
# lst.mat.grn.cns.comparison <- vector(mode = "list", length = 2)
# 
# print("characteristics of equal size comparison network")
# idx.equal_size <- max(which(df.rate_density$n.pred <= lst.eval$n.pred))
# print(df.rate_density[idx.equal_size,])
# 
# 
# th.grn.comparison <- df.rate_density$v.th.grn[idx.equal_size]
# lst.mat.grn.cns.comparison[[1]] <- mat.gene_regulatory_network
# lst.mat.grn.cns.comparison[[1]][lst.mat.grn.cns.comparison[[1]] < th.grn.comparison] <- 0
# lst.mat.grn.cns.comparison[[1]][lst.mat.grn.cns.comparison[[1]] >= th.grn.comparison] <- 1
# 
# mat.grn.cns.comparison=lst.mat.grn.cns.comparison[[1]]
# 
# # general overlap and network properties
# compute_NetworkOverlap(mat.p.grace_ensemble, mat.grn.cns.comparison, trace = TRUE, test = "hyper")
# print(paste("GRACE: # regulators:", length(which(rowSums(mat.p.grace_ensemble) > 0)), " # targets:", length(which(colSums(mat.p.grace_ensemble) > 0)) ))
# print(paste("COMPARISON: # regulators:", length(which(rowSums(mat.grn.cns.comparison) > 0)), " # targets:", length(which(colSums(mat.grn.cns.comparison) > 0)) ))
# 
# # coregulation enrichment
# hitInSample <- lst.eval$n.coreg.tg_pairs
# sampleSize <-  length(which(colSums(mat.p.grace_ensemble) > 0)) *  length(which(colSums(mat.p.grace_ensemble) > 0))
# hitInPop <- df.rate_density$n.coreg.tg_pairs[idx.equal_size]
# popSize <- length(which(colSums(mat.grn.cns.comparison) > 0)) * length(which(colSums(mat.grn.cns.comparison) > 0))
# (hitInSample/sampleSize)/(hitInPop/popSize)
# 
# tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
# fisher.test(tab)$p.value
# 
# ### validation double annotation exclusion (make automatic)
# # 
# tf.gs <- intersect(rownames(mat.p.grace_ensemble), rownames(lst.benchmarks$mat.regulatory_evidence))
# tf.gs <- intersect(tf.gs, rownames(mat.AtRegNet))
# tg.gs <- intersect(colnames(mat.p.grace_ensemble), colnames(lst.benchmarks$mat.regulatory_evidence))
# tg.gs <- intersect(tg.gs, colnames(mat.AtRegNet))
# 
# sum(mat.p.grace_ensemble[tf.gs,tg.gs] * mat.AtRegNet[tf.gs,tg.gs])
# 
# sum(mat.p.grace_ensemble[tf.gs,tg.gs] * lst.benchmarks$mat.regulatory_evidence[tf.gs,tg.gs])
# 
# test <- (mat.AtRegNet[tf.gs,tg.gs] * lst.benchmarks$mat.regulatory_evidence[tf.gs,tg.gs]) # overlap of AGRIS and ATRM
# 
# sum(test * mat.p.grace_ensemble[tf.gs,tg.gs]) # present in both datasets
# 
# ###
# 
# mat.grn.MR <- jaccard(t(as.matrix(mat.p.grace_ensemble)))
# mat.grn.MR <- as.matrix(mat.grn.MR); diag(mat.grn.MR) <- 0;
# rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(mat.p.grace_ensemble) 
# mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
# 
# mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
# mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
# 
# tg.gs <- intersect(colnames(mat.grn.MR), colnames(lst.benchmarks$mat.cofunctional_evidence))
# tg.gs <- intersect(tg.gs, colnames(mat.aracyc.tgs))
# 
# sum(mat.grn.MR[tg.gs,tg.gs] * mat.aracyc.tgs[tg.gs,tg.gs])
# sum(mat.grn.MR[tg.gs,tg.gs] * lst.benchmarks$mat.cofunctional_evidence[tg.gs,tg.gs])
# 
# test <- (mat.aracyc.tgs[tg.gs,tg.gs] * lst.benchmarks$mat.cofunctional_evidence[tg.gs,tg.gs])
# 
# sum(test * mat.grn.MR[tg.gs,tg.gs]) # present in both datasets (Aracyc)
# 
# ### comparison network
# 
# tf.gs <- intersect(rownames(mat.grn.cns.comparison), rownames(lst.benchmarks$mat.regulatory_evidence))
# tf.gs <- intersect(tf.gs, rownames(mat.AtRegNet))
# tg.gs <- intersect(colnames(mat.grn.cns.comparison), colnames(lst.benchmarks$mat.regulatory_evidence))
# tg.gs <- intersect(tg.gs, colnames(mat.AtRegNet))
# 
# sum(mat.grn.cns.comparison[tf.gs,tg.gs] * mat.AtRegNet[tf.gs,tg.gs])
# sum(mat.grn.cns.comparison[tf.gs,tg.gs] * lst.benchmarks$mat.regulatory_evidence[tf.gs,tg.gs])
# 
# test <- (mat.AtRegNet[tf.gs,tg.gs] * lst.benchmarks$mat.regulatory_evidence[tf.gs,tg.gs])
# sum(test)
# sum(test * mat.grn.cns.comparison[tf.gs,tg.gs]) # present in both datasets (agris vs atrm)
# 
# 
# 
# 
# mat.grn.MR <- jaccard(t(as.matrix(lst.mat.grn.cns.comparison[[1]])))
# mat.grn.MR <- as.matrix(mat.grn.MR); diag(mat.grn.MR) <- 0;
# rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(lst.mat.grn.cns.comparison[[1]]) 
# mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
# 
# mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
# mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
# 
# tg.gs <- intersect(colnames(mat.grn.MR), colnames(mat.aracyc.tgs))
# tg.gs <- intersect(colnames(mat.grn.MR), colnames(lst.benchmarks$mat.cofunctional_evidence))
# tg.gs <- intersect(tg.gs, colnames(mat.aracyc.tgs))
# 
# sum(mat.grn.MR[tg.gs,tg.gs] * mat.aracyc.tgs[tg.gs,tg.gs])
# sum(mat.grn.MR[tg.gs,tg.gs] * lst.benchmarks$mat.cofunctional_evidence[tg.gs,tg.gs])
# 
# test <- (mat.aracyc.tgs[tg.gs,tg.gs] * lst.benchmarks$mat.cofunctional_evidence[tg.gs,tg.gs])
# 
# sum(test * mat.grn.MR[tg.gs,tg.gs]) # present in both datasets
# 
# ##
# 
# print("characteristics of best fscore comparison network")
# idx.best_fscore <- max(which(df.rate_density$fscore_beta_01 == max(df.rate_density$fscore_beta_01)))
# print(df.rate_density[idx.best_fscore,])
# 
# th.grn.comparison <- df.rate_density$v.th.grn[idx.best_fscore]
# lst.mat.grn.cns.comparison[[2]] <- mat.gene_regulatory_network
# lst.mat.grn.cns.comparison[[2]][lst.mat.grn.cns.comparison[[2]] < th.grn.comparison] <- 0
# lst.mat.grn.cns.comparison[[2]][lst.mat.grn.cns.comparison[[2]] >= th.grn.comparison] <- 1
# 
# ## double check
# tf.gs <- intersect(rownames(lst.mat.grn.cns.comparison[[2]]), rownames(lst.benchmarks$mat.regulatory_evidence))
# tf.gs <- intersect(tf.gs, rownames(mat.AtRegNet))
# tg.gs <- intersect(colnames(lst.mat.grn.cns.comparison[[2]]), colnames(lst.benchmarks$mat.regulatory_evidence))
# tg.gs <- intersect(tg.gs, colnames(mat.AtRegNet))
# 
# sum(lst.mat.grn.cns.comparison[[2]][tf.gs,tg.gs] * mat.AtRegNet[tf.gs,tg.gs])
# sum(lst.mat.grn.cns.comparison[[2]][tf.gs,tg.gs] * lst.benchmarks$mat.regulatory_evidence[tf.gs,tg.gs])
# 
# test <- (mat.AtRegNet[tf.gs,tg.gs] * lst.benchmarks$mat.regulatory_evidence[tf.gs,tg.gs])
# sum(test)
# sum(test * lst.mat.grn.cns.comparison[[2]][tf.gs,tg.gs]) # present in both datasets
# 
# 
# mat.grn.MR <- jaccard(t(as.matrix(lst.mat.grn.cns.comparison[[2]])))
# mat.grn.MR <- as.matrix(mat.grn.MR); diag(mat.grn.MR) <- 0;
# rownames(mat.grn.MR) <- colnames(mat.grn.MR) <- colnames(lst.mat.grn.cns.comparison[[2]]) 
# mat.grn.MR <- as(mat.grn.MR, "CsparseMatrix")
# 
# mat.grn.MR@x[mat.grn.MR@x < 0.5] <- 0
# mat.grn.MR@x[mat.grn.MR@x >= 0.5] <- 1
# 
# tg.gs <- intersect(colnames(mat.grn.MR), colnames(mat.aracyc.tgs))
# tg.gs <- intersect(colnames(mat.grn.MR), colnames(lst.benchmarks$mat.cofunctional_evidence))
# tg.gs <- intersect(tg.gs, colnames(mat.aracyc.tgs))
# 
# sum(mat.grn.MR[tg.gs,tg.gs] * mat.aracyc.tgs[tg.gs,tg.gs])
# sum(mat.grn.MR[tg.gs,tg.gs] * lst.benchmarks$mat.cofunctional_evidence[tg.gs,tg.gs])
# 
# test <- (mat.aracyc.tgs[tg.gs,tg.gs] * lst.benchmarks$mat.cofunctional_evidence[tg.gs,tg.gs])
# 
# sum(test * mat.grn.MR[tg.gs,tg.gs]) # present in both datasets
# 
# 
# # for paper only - not part of the software
# barplot <- (evaluate_GRACE(lst.eval_grace=lst.eval, v.comparison_best_fscore = df.rate_density[idx.best_fscore,], v.comparison_equal_size = df.rate_density[idx.equal_size,], v.background = df.rate_density[nrow(df.rate_density),],
#                            v.network_sizes = c("2192","792","792"),
#                            v.compensation.agris = c(3,2,7), 
#                            v.compensation.subacon = c(165,18,32),
#                            v.compensation.aracyc = c(33,0,7))) 
# 
# 
# #cairo_pdf("example.pdf", width = 7.2, height = 7, family= "Helvetica")#"DejaVu Sans")
# plot(barplot)
# #dev.off()
# 
# 

