


# optimiation function - grid search
gridSearchOptimization <- function(lst.gs.grn_bp.train=lst.gs.grn_bp.train, gamma.gridSearch, lambda.gridSearch, th.percentage.hyperparameters = 0.99, ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization, max.coreg.bp = max.coreg.bp){
  
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
      
      mat.p.crf <- lst.mat.p.crf$mat.p.crf # get if equal to zero ... , could be based on bootstrap sample, no value in there... 
      
      # use the fast model
      lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                              lst.benchmarks=lst.gs.grn_bp.train, 
                                                              ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization, 
                                                              max.coreg.bp = max.coreg.bp)
                                                        
      #       lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
      #                                                   lst.benchmarks=lst.gs.grn_bp.train,
      #                                                   th.min.coreg.tfs = th.min.coreg.tfs, 
      #                                                   b.jaccard = b.jaccard,
      #                                                   beta = beta)
            
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
    
    # use the fast, simpler function
    lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                            lst.benchmarks=lst.gs.grn_bp.train)
    
    #     lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
    #                                                 lst.benchmarks=lst.gs.grn_bp.train,
    #                                                 th.min.coreg.tfs = th.min.coreg.tfs, 
    #                                                 b.jaccard = b.jaccard,
    #                                                 beta = beta)
    #     
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


# eval 
evaluate_mode_parameter <- function(lst.grace_support=lst.grace_support, 
                                    lst.modules=lst.modules, 
                                    truereg=truereg, 
                                    mat.cofunctional_evidence=mat.cofunctional_evidence, 
                                    
                                    gamma.gridSearch=gamma.gridSearch,
                                    lambda.gridSearch=lambda.gridSearch,
                                    th.percentage.hyperparameters=th.percentage.hyperparameters,
                                    
                                    ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization,
                                    max.coreg.bp = max.coreg.bp,
                                    
                                    n.sample_size=n.sample_size, 
                                    max.call = 100, 
                                    v.seed = 1234, 
                                    beta = 1.0){
  
  source("GRACE_optimization_routines.R")
  source("GRACE.R")
  source("GRACE_helperfunctions.R")

  #### Bootstrapping (based on regulators based meta module sampling (0.632))
  
  set.seed(v.seed)
  
  tfs.mrf <- unique(lst.grace_support$df.grn$TF) 
  tgs.mrf <- unique(lst.grace_support$df.grn$TG)
  
  lst.bootstrap_sample <- generate_bootstrapping_sample(lst.grace_support=lst.grace_support,lst.modules=lst.modules, truereg=truereg, mat.cofunctional_evidence=mat.cofunctional_evidence, n.sample_size = n.sample_size)
  lst.grace_support.bootstrap_sample <- lst.bootstrap_sample$lst.grace_support.bootstrap_sample
  lst.modules.bootstrap_sample <- lst.bootstrap_sample$lst.modules.bootstrap_sample
  
  lst.gs.grn_bp.train <- list(lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.GS.bootstrap_sample, lst.bootstrap_sample$lst.benchmark.bootstrap_sample$mat.cofunctional_evidence.bootstrap_sample)    #list(mat.GS.train, mat.BP.coreg.train)
  names(lst.gs.grn_bp.train) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
  
  ## Step 1: Grid Search (Coarse) optimization to identify finer range for Simulated annealing optimization ##
  strt<-Sys.time() 
  #lst.hyperparameter <- gridSearchOptimization(lst.gs.grn_bp.train=lst.gs.grn_bp.train,gamma.gridSearch, lambda.gridSearch, th.percentage.hyperparameters, ratio.gs.grn_vs_bp.normalization = ratio.gs.grn_vs_bp.normalization, max.coreg.bp = max.coreg.bp)
  
  # add begin grid search
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
          
          mat.p.crf <- lst.mat.p.crf$mat.p.crf # get if equal to zero ... , could be based on bootstrap sample, no value in there... 
          
          if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
            
            lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                                    lst.benchmarks=lst.gs.grn_bp.train)
            
            
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
    
      #return(list(w=w,lower_limits=lower_limits,upper_limits=upper_limits))

      w.intials     <- w
      lower_limits <- lower_limits
      upper_limits <- upper_limits
      
        
      #   w.intials     <- lst.hyperparameter$w
      #   lower_limits <- lst.hyperparameter$lower_limits
      #   upper_limits <- lst.hyperparameter$upper_limits
      #   
        
  print(Sys.time()-strt)
  
#   #     w     <- c(gamma.upper, 1.5) #0.5) # start at upper end - therefore, the one with smallest 
#   #     lower <- c(gamma.lower, 2.0) #0.01) 
#   #     upper <- c(gamma.upper, 2.5) # 1.0)
#   
  ## Step2: Simulated annealing ##

  
  
  strt<-Sys.time() 
  
  # begin 
    evaluate_SA <- function(w) {
      
      v.alpha = v.gamma = w[1]; v.lambda= w[2]; 
      if(v.alpha >= v.gamma){
        lst.mat.p.crf <- reevaluate_links_with_MarkovRandomFields (lst.grace_support=lst.grace_support.bootstrap_sample, lst.modules=lst.modules.bootstrap_sample,
                                                                   v.alpha = v.alpha, v.gamma = v.gamma, v.lambda = v.lambda, b.prim=TRUE, tfs = tfs.mrf, tgs = tgs.mrf) 
        
        mat.p.crf <- lst.mat.p.crf$mat.p.crf
        
        if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
          
        
          # use the fast, simpler function
          lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                                  lst.benchmarks=lst.gs.grn_bp.train)
          
          #     lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
          #                                                 lst.benchmarks=lst.gs.grn_bp.train,
          #                                                 th.min.coreg.tfs = th.min.coreg.tfs, 
          #                                                 b.jaccard = b.jaccard,
          #                                                 beta = beta)
          
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
  ## end define SA
  
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
  
  if(!is.null(mat.p.crf) & !is.null(dim(mat.p.crf)) & sum(mat.p.crf) > 0){
    
    lst.gs.grn_bp <- list(truereg, mat.cofunctional_evidence)
    names(lst.gs.grn_bp) <- c("mat.regulatory_evidence", "mat.cofunctional_evidence")
    
    lst.eval <- compute_model_assessment_parallel_optimized(mat.p.grace_ensemble=mat.p.crf,
                                                            lst.benchmarks=lst.gs.grn_bp)
    
    #     lst.eval <- compute_fscore_model_assessment(mat.p.grace_ensemble=mat.p.crf,
    #                                                 lst.benchmarks=lst.gs.grn_bp,  
    #                                                 th.min.coreg.tfs = th.min.coreg.tfs, 
    #                                                 b.jaccard = b.jaccard,
    #                                                 beta = beta)
    
    
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
  }else{
    lst.grace_models <- list(weights=c(-1, -1, -1), fscore=-100, n.pred=-1,n.gs.grn=-1,n.coreg.bp.tgs=-1)  
  }
  
  return(lst.grace_models)
}

