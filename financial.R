
library(glmnet)
library(foreach)
library(doParallel)

source("TL-FLR-FunFiles.R")

datafile <- c('May_to_Jun', 'Jun_to_Jul', 'Jul_to_Aug', 'Aug_to_Sept')

ncores <- 25
mc <- 500

split_ratio <- 0.8
# candidate tuning parameters
lamseq <- c(exp(seq(log(5), log(1e-3), length.out = 10)))
mseq <- 1:8
params <- list(omega=c(50, 100, 250, 1250), niter=100, tol=1e-6) 

for(idx_data in 1:4){
  load(paste0('Financial/', datafile[idx_data], '.RData'))
  L <- length(sectors_data_list)
  sec_names <- names(sectors_data_list)

  for(i in 1:L){
    # target
    name_target <- sec_names[i]
    x <- sectors_data_list[[i]]$X
    y <- sectors_data_list[[i]]$Y
    n <- nrow(x)
    
    # source
    idx <- setdiff((1:L), i)
    x.S <- vector(mode='list', length=L-1)
    y.S <- vector(mode='list', length=L-1)
    for(j in 1:(L-1)){
      x.S[[j]] <- sectors_data_list[[idx[j]]]$X
      y.S[[j]] <- sectors_data_list[[idx[j]]]$Y
    }
    n.S <- sapply(y.S, length)
    
    nt <- ncol(x)
    tseq <- seq(0, 1, length.out=nt)
    t.by <- tseq[2]-tseq[1]
    
    cl <- makeCluster(ncores)  
    registerDoParallel(cl)
    rec <- foreach(run = 1:mc, .packages = c("glmnet")) %dopar% {
      ntrain <- ceiling(split_ratio*n); ntest <- n - ntrain
      idx_train <- sample(n, ntrain)
      xtrain <- x[idx_train, ]; xtest <- x[-idx_train, ]
      ytrain <- y[idx_train]; ytest <- y[-idx_train]
      
      ###############################
      # using only target data
      # the tuning parameter is chosen by CV
      ###############################
      cv_result <- cv.FLR(xtrain, ytrain, NULL, NULL, tseq, mseq, k=5)
      m.cv_flr <- cv_result$m.opt
      result <- FLR(xtrain, ytrain, NULL, NULL, tseq, m.cv_flr)
      slope.est <- result$slope.est
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
      preerror.cv_flr <- mean((ytest - a.est - xtest%*%slope.est*t.by)^2)
      
      ###############################
      # naive transfer learning: all source samples
      # the tuning parameter is chosen by 5-fold CV
      ###############################
      result_step1 <- tlFLR_step1(xtrain, ytrain, x.S, y.S, tseq, max(mseq))
      cv_result <- cv.tlFLR(xtrain, ytrain, result_step1$w.est, result_step1$evec.est, tseq, mseq, lamseq, k=5)
      m.naive_cv_tl_flr <- cv_result$m.opt
      lam.naive_cv_tl_flr <- cv_result$lam.opt
      # result_debias <- tlFLR_step2(xtrain, ytrain, tseq, m.naive_cv_tl_flr, lam.naive_cv_tl_flr, result_step1$w.est, result_step1$evec.est)
      # slope.est <- result_debias$slope.est
      if (m.naive_cv_tl_flr==1) {
        result_debias <- tlFLR_step2(xtrain, ytrain, tseq, m.naive_cv_tl_flr, 0, result_step1$w.est, result_step1$evec.est)
        slope.est <- result_debias$slope.est
      }else {
        result_debias <- tlFLR_step2(xtrain, ytrain, tseq, m.naive_cv_tl_flr, lamseq, result_step1$w.est, result_step1$evec.est)
        idx_lam <- which(lamseq==lam.naive_cv_tl_flr)
        slope.est <- result_debias$slope.est[,idx_lam]
      }
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
      preerror.naive_cv_tl_flr <- mean((ytest - a.est - xtest%*%slope.est*t.by)^2)
      
      ###############################
      # aggregation transfer learning; contrast: g function
      ###############################
      result <- tlFLR_Agg(xtrain, ytrain, x.S, y.S, tseq, mseq, lamseq, params, num_splits=5, split_ratio=0.6, k=5)
      slope.est.qagg <- result$slope.est.qagg
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est.qagg*t.by)
      preerror.cv_qagg_tl_flr <- mean((ytest - a.est - xtest%*%slope.est.qagg*t.by)^2)
      
      slope.est.sagg <- result$slope.est.sagg
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est.sagg*t.by)
      preerror.cv_sagg_tl_flr <- mean((ytest - a.est - xtest%*%slope.est.sagg*t.by)^2)
      
      ###############################
      # aggregation pooled estimate (without de-bias); contrast: g function
      ###############################
      result <- FLR_Agg(xtrain, ytrain, x.S, y.S, tseq, mseq, params, num_splits=5, split_ratio=0.6, k=5)
      slope.est.qagg <- result$slope.est.qagg
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est.qagg*t.by)
      preerror.cv_qagg_pool_flr <- mean((ytest - a.est - xtest%*%slope.est.qagg*t.by)^2)
      
      slope.est.sagg <- result$slope.est.sagg
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est.sagg*t.by)
      preerror.cv_sagg_pool_flr <- mean((ytest - a.est - xtest%*%slope.est.sagg*t.by)^2)
      
      return(list(preerror.cv_flr=preerror.cv_flr, preerror.naive_cv_tl_flr=preerror.naive_cv_tl_flr,
                  preerror.cv_qagg_tl_flr=preerror.cv_qagg_tl_flr, preerror.cv_sagg_tl_flr=preerror.cv_sagg_tl_flr,
                  preerror.cv_qagg_pool_flr=preerror.cv_qagg_pool_flr, preerror.cv_sagg_pool_flr=preerror.cv_sagg_pool_flr))
    }
    stopCluster(cl)
    save(rec, file=paste0(datafile[idx_data], '_Error_', name_target, '.Rdata'))
  }
}


