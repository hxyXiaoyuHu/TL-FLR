
library(glmnet)
library(foreach)
library(doParallel)

source("TL-FLR-FunFiles.R")

# model parameters
n <- 150
alpha <- 1.5
beta <- 2
sigma.e <- 0.5
nt <- 100
tseq <- seq(0.01, 1, length.out=nt)
t.by <- tseq[2]-tseq[1]
L <- 20; s <- 50
n.S <- rep(100, L)
Kseq <- seq(0, 20, by=2)
Klen <- length(Kseq)

ncores <- 20
mc <- 500
sourceBasisType <- 'haar'
# candidate tuning parameters
lamseq <- c(exp(seq(log(1), log(1e-3), length.out = 10)))
mseq <- seq(1, 10, length.out=10)
params <- list(omega=c(0.2, 0.5, 2, 10), niter=100, tol=1e-6)

hseq <- c(2, 20, 200)
hlen <- length(hseq)
for(idx_h in 1:3){
  h <- hseq[idx_h]
  D <- D.S <- diag((1:s)^(-alpha/2))
    
  for (idx_K in 1:Klen) {
    K <- Kseq[idx_K]
    print(paste0("the number of informative sets: ", K))
    
    n.S <- rep(100, L)
    b <- 4*(1:s)^(-beta)*(-1)^((1:s)+1)
    b.S <- matrix(0, s, L)
   
    R <- sample(c(-1,1), size=s, replace=TRUE)
    for (l in 1:L) {
      if (l <= K) { # informative
        b.S[,l] <- b - h/s*R
      } else { # non-informative
        b.S[,l] <- b - 40*R
      }
    }        
      
    cl <- makeCluster(ncores)  
    registerDoParallel(cl)
    rec <- foreach(run = 1:mc, .packages = c("glmnet")) %dopar% {
      dat <- genData(n, n.S, D, D.S, tseq, b, b.S, K=K, sigma.e=sigma.e, s=s, sourceBasisType=sourceBasisType)
      x <- dat$x; y <- dat$y
      x.S <- dat$x.S; y.S <- dat$y.S
      slope.val <- dat$slope.val
 
      ###############################
      # using only target data
      # the tuning parameter is chosen by CV
      ###############################
      cv_result <- cv.FLR(x, y, NULL, NULL, tseq, mseq, k=5)
      m.cv_flr <- cv_result$m.opt
      result <- FLR(x, y, NULL, NULL, tseq, m.cv_flr)
      slope.est <- result$slope.est
      error.cv_flr <- sum((slope.est - slope.val)^2 * t.by) 
      
      ###############################
      # naive transfer learning: all source samples
      # the tuning parameter is chosen by 5-fold CV
      ###############################
      result_step1 <- tlFLR_step1(x, y, x.S, y.S, tseq, max(mseq))
      cv_result <- cv.tlFLR(x, y, result_step1$w.est, result_step1$evec.est, tseq, mseq, lamseq, k=5)
      m.naive_cv_tl_flr <- cv_result$m.opt
      lam.naive_cv_tl_flr <- cv_result$lam.opt
      # result_debias <- tlFLR_step2(x, y, tseq, m.naive_cv_tl_flr, lam.naive_cv_tl_flr, result_step1$w.est, result_step1$evec.est)
      # slope.est <- result_debias$slope.est
      # error.naive_cv_tl_flr <- sum((slope.est - slope.val)^2 * t.by)
      if (m.naive_cv_tl_flr==1) {
        result_debias <- tlFLR_step2(x, y, tseq, m.naive_cv_tl_flr, 0, result_step1$w.est, result_step1$evec.est)
        slope.est <- result_debias$slope.est
      }else {
        result_debias <- tlFLR_step2(x, y, tseq, m.naive_cv_tl_flr, lamseq, result_step1$w.est, result_step1$evec.est)
        idx_lam <- which(lamseq==lam.naive_cv_tl_flr)
        slope.est <- result_debias$slope.est[,idx_lam]
      }
      error.naive_cv_tl_flr <- sum((slope.est - slope.val)^2 * t.by)
      
      ###############################
      # oracle transfer learning: known informative sets
      # the tuning parameter is chosen by 5-fold CV
      ###############################
      if(K==0) {
        m.cv_tl_flr <- m.cv_flr
        lam.cv_tl_flr <- 0
        error.cv_tl_flr <- error.cv_flr
      }else {
        result_step1 <- tlFLR_step1(x, y, x.S[1:K], y.S[1:K], tseq, max(mseq))
        cv_result <- cv.tlFLR(x, y, result_step1$w.est, result_step1$evec.est, tseq, mseq, lamseq, k=5)
        m.cv_tl_flr <- cv_result$m.opt
        lam.cv_tl_flr <- cv_result$lam.opt
        # result_debias <- tlFLR_step2(x, y, tseq, m.cv_tl_flr, lam.cv_tl_flr, result_step1$w.est, result_step1$evec.est)        
        # slope.est <- result_debias$slope.est
        # error.cv_tl_flr <- sum((slope.est - slope.val)^2 * t.by)
        if (m.cv_tl_flr==1) {
          result_debias <- tlFLR_step2(x, y, tseq, m.cv_tl_flr, 0, result_step1$w.est, result_step1$evec.est)
          slope.est <- result_debias$slope.est
        } else {
          result_debias <- tlFLR_step2(x, y, tseq, m.cv_tl_flr, lamseq, result_step1$w.est, result_step1$evec.est)
          idx_lam <- which(lamseq == lam.cv_tl_flr)
          slope.est <- result_debias$slope.est[,idx_lam]
        }
        error.cv_tl_flr <- sum((slope.est - slope.val)^2 * t.by)
      }
      
      ###############################
      # aggregation transfer learning; contrast: g function
      ###############################
      result <- tlFLR_Agg(x, y, x.S, y.S, tseq, mseq, lamseq, params, num_splits=5, split_ratio=0.6, k=5)
      slope.est.qagg <- result$slope.est.qagg
      error.cv_qagg_tl_flr <- sum((slope.est.qagg - slope.val)^2 * t.by)
      
      slope.est.sagg <- result$slope.est.sagg
      error.cv_sagg_tl_flr <- sum((slope.est.sagg - slope.val)^2 * t.by)
      
      ###############################
      # aggregation pooled estimate (without de-bias); contrast: g function
      ###############################
      result <- FLR_Agg(x, y, x.S, y.S, tseq, mseq, params, num_splits=5, split_ratio=0.6, k=5)
      slope.est.qagg <- result$slope.est.qagg
      error.cv_qagg_pool_flr <- sum((slope.est.qagg - slope.val)^2 * t.by)
      
      slope.est.sagg <- result$slope.est.sagg
      error.cv_sagg_pool_flr <- sum((slope.est.sagg - slope.val)^2 * t.by)
      
      return(list(error.cv_flr = error.cv_flr, m.cv_flr=m.cv_flr, 
                  error.cv_tl_flr=error.cv_tl_flr,
                  m.cv_tl_flr=m.cv_tl_flr, lam.cv_tl_flr=lam.cv_tl_flr,
                  error.cv_qagg_tl_flr=error.cv_qagg_tl_flr, error.cv_sagg_tl_flr=error.cv_sagg_tl_flr,
                  error.cv_qagg_pool_flr=error.cv_qagg_pool_flr, error.cv_sagg_pool_flr=error.cv_sagg_pool_flr,
                  error.naive_cv_tl_flr=error.naive_cv_tl_flr,
                  m.naive_cv_tl_flr=m.naive_cv_tl_flr, lam.naive_cv_tl_flr=lam.naive_cv_tl_flr))
    }
    stopCluster(cl)
    save(rec, file=paste0("rec_", sourceBasisType, "_h_", h, "_K_", K, ".Rdata"))
  }
}
