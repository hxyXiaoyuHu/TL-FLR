
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
mc <- 1000
sourceBasisType <- 'fourier'
# candidate tuning parameters
lamseq <- c(exp(seq(log(1), log(1e-3), length.out = 10)))
mseq <- seq(1, 10, length.out=10)

hseq <- c(2, 20, 200, 2000)
hlen <- length(hseq)
ss_seq <- c(1, 5, 20, 50)
for(idx_ss in 1:4){
  ss <- ss_seq[idx_ss]
  for(idx_h in 1:4){
    h <- hseq[idx_h]
    D <- D.S <- diag((1:s)^(-alpha/2))
    
    for (idx_K in 2:Klen) {
      K <- Kseq[idx_K]
      L <- K # known informative sets (h-level)
      
      n.S <- rep(100, L)
      b <- 4*(1:s)^(-beta)*(-1)^((1:s)+1)
      b.S <- matrix(0, s, L)
     
      R <- sample(c(-1,1), size=s, replace=TRUE)
      for (l in 1:L) {
        b.S[,l] <- b     
        b.S[1:ss,l] <- b[1:ss] - h/ss*R[1:ss]
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
        # naively using both target and source data
        # the tuning parameter is chosen by 5-fold CV
        ###############################
        cv_result <- cv.FLR(x, y, x.S, y.S, tseq, mseq, k=5)
        m.cv_pool_flr <- cv_result$m.opt
        result <- FLR(x, y, x.S, y.S, tseq, m.cv_pool_flr)
        slope.est <- result$slope.est
        error.cv_pool_flr <- sum((slope.est - slope.val)^2 * t.by)
        
        ###############################
        # oracle transfer learning: known informative sets
        # the tuning parameter is chosen by 5-fold CV
        ###############################
        result_step1 <- tlFLR_step1(x, y, x.S, y.S, tseq, max(mseq))
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
          slope.est <- result_debias$slope.est[, idx_lam]
        }
        error.cv_tl_flr <- sum((slope.est - slope.val)^2 * t.by)
        
        return(list(error.cv_flr = error.cv_flr, error.cv_pool_flr=error.cv_pool_flr, 
                    error.cv_tl_flr=error.cv_tl_flr, 
                    m.cv_flr=m.cv_flr, m.cv_pool_flr=m.cv_pool_flr,
                    m.cv_tl_flr=m.cv_tl_flr, lam.cv_tl_flr=lam.cv_tl_flr))
      }
      stopCluster(cl)
      save(rec, file=paste0("rec_ndiffcoef_", ss, "_h_", h, "_K_", K, ".Rdata"))
    }
  }
}
