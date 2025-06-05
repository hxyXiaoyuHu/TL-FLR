
# @param n, sample size of target sample
# @param n.S, sample sizes of source samples, length L
# @param D/D.S, square root of diagonal covariance matrix, s*s
# @param tseq, sampling points, length nt (regular design)
# @param b/b.S, coefficients of target/source data, length s/matrix s*L
# @param sigma.e, standard deviation of random noise
# @param s, number of basis functions; 
# @param L, number of source samples; K, number of informative source samples

library(glmnet)

fourier.basis <- function(s, tseq) {
  nt <- length(tseq)
  basis.mat <- matrix(NA, nrow=nt, ncol=s)
  for (j in 1:s) {
    basis.mat[,j] <- sqrt(2)*cos(j*pi*tseq)
  }
  return(basis.mat)
}

# Define the Haar function
haar_function <- function(k, t) {
  if (k == 0) return(1)
  
  j <- floor(log2(k))
  m <- k - 2^j
  
  scale_factor <- 2^(j / 2)
  left <- m / 2^j
  middle <- (m + 0.5) / 2^j
  right <- (m + 1) / 2^j
  
  if (t >= left && t < middle) {
    return(scale_factor)
  } else if (t >= middle && t <= right) {
    return(-scale_factor)
  } else {
    return(0)
  }
}

# Generate Haar function values
generate_haar_values <- function(k, x_values) {
  sapply(x_values, function(x) haar_function(k, x))
}

haar.basis <- function(s, tseq){
  nt <- length(tseq)
  basis.mat <- matrix(NA, nt, s)
  for (k in 1:s) {
    basis.mat[,k] <- generate_haar_values(k, tseq)
  }
  return(basis.mat)
}


genData <- function(n, n.S, D, D.S, tseq, b, b.S, K, sigma.e=0.5, s=50, sourceBasisType='fourier') {
  
  # target data
  basis.mat <- fourier.basis(s,tseq)
  z <- matrix(runif(n*s, -sqrt(3), sqrt(3)), n, s) # uniform distribution
  # z <- sqrt(3/5) * matrix(rt(n*s, df=5), n, s) # t(5)
  score <- z%*%D
  x <- tcrossprod(score, basis.mat)
  slope.val <- basis.mat%*%b
  y <- score%*%b + rnorm(n, sd=sigma.e)
  
  # source data
  L <- length(n.S)
  if (L==0) { # there is no source data
    x.S <- NULL 
    y.S <- NULL
  } else {
    if (sourceBasisType=='fourier'){
      x.S <- vector(mode="list", length=L)
      y.S <- vector(mode="list", length=L)
      for (l in 1:L) {
        nl <- n.S[l]
        z <- matrix(rnorm(nl*s, 0, 1), nl, s) # normal distribution
        score <- z%*%D.S
        x.S[[l]] <- tcrossprod(score, basis.mat)
        y.S[[l]] <- score%*%b.S[,l] + rnorm(nl, sd=sigma.e)
      }
    } else if (sourceBasisType=='haar') {
      t.by <- tseq[2] - tseq[1]
      basis.mat.S <- haar.basis(s, tseq)
      x.S <- vector(mode="list", length=L)
      y.S <- vector(mode="list", length=L)
      for (l in 1:L) {
        nl <- n.S[l]
        z <- matrix(rnorm(nl*s, 0, 1), nl, s) # normal distribution
        score <- z%*%D.S
        x.S[[l]] <- tcrossprod(score, basis.mat.S)
        y.S[[l]] <- x.S[[l]] %*% basis.mat %*% b.S[,l] * t.by + rnorm(nl, sd=sigma.e)
      }
    } else if (sourceBasisType=='fourier-haar'){
      t.by <- tseq[2] - tseq[1]
      basis.mat.S <- haar.basis(s, tseq)
      x.S <- vector(mode="list", length=L)
      y.S <- vector(mode="list", length=L)
      for (l in 1:L) {
        nl <- n.S[l]
        z <- matrix(rnorm(nl*s, 0, 1), nl, s) # normal distribution
        score <- z%*%D.S
        if(l <= K){
          x.S[[l]] <- tcrossprod(score, basis.mat) # Fourier basis
          y.S[[l]] <- score%*%b.S[,l] + rnorm(nl, sd=sigma.e)
        }else {
          x.S[[l]] <- tcrossprod(score, basis.mat.S) # Haar basis
          y.S[[l]] <- x.S[[l]] %*% basis.mat %*% b.S[,l] * t.by + rnorm(nl, sd=sigma.e)
        }
      }
    } 
  }
  
  return(list(x=x, y=y, x.S=x.S, y.S=y.S, slope.val=slope.val))
}

estFPCA <- function(x, tseq, m=NULL, FVE=0.95) {
  # @param x, list
  
  L <- length(x)
  nt <- length(tseq); t.by <- tseq[2] - tseq[1]
  n.S <- sapply(x, nrow); N <- sum(n.S)
  cov.S <- 0
  for (l in 1:L) {
    xl <- x[[l]]; nl <- n.S[l]
    cov.S <- nl/N * cov(xl) + cov.S # unbiased estimate
  }
  eig.est <- eigen(cov.S)
  eval.est <- eig.est$values*t.by
  eval.est <- eval.est[eval.est>0]
  if(is.null(m)){
    n.pc <- min(which((cumsum(eval.est)/sum(eval.est))>=FVE))
  }else{
    if(length(eval.est)>=m){
      n.pc <- m
    }else{
      n.pc <- length(eval.est)
    } 
  }
  eval.est <- eval.est[1:n.pc]
  evec.est <- eig.est$vectors[,1:n.pc,drop=F]/sqrt(t.by)
  
  return(list(evec.est=evec.est, eval.est=eval.est))
}

listConcatenate <- function(x, y, x.S, y.S){
  L <- length(x.S)
  x.all <- vector(mode='list', length = L+1)
  y.all <- vector(mode='list', length = L+1)
  for (l in 1:(L+1)) {
    if (l==1) {
      x.all[[l]] <- x
      y.all[[l]] <- y
    } else {
      x.all[[l]] <- x.S[[l-1]]
      y.all[[l]] <- y.S[[l-1]]
    }
  }
  return(list(x=x.all, y=y.all))
}

tlFLR_step1 <- function(x, y, x.S, y.S, tseq, m, targetInc=TRUE) {
  # transfer learning algorithm: step 1
  # biased estimate using source datasets
  # if m is a vector, then set the truncation number to max(m)
  
  if(targetInc){
    tmp <- listConcatenate(x, y, x.S, y.S)
    x.S <- tmp$x
    y.S <- tmp$y
  }
  if(is.null(m)){
    result.FPCA <- estFPCA(x.S, tseq)
  }else{
    result.FPCA <- estFPCA(x.S, tseq, m=max(m))
  }
  evec.est <- result.FPCA$evec.est
  eval.est <- result.FPCA$eval.est
  n.pc <- ncol(evec.est)
  
  L <- length(y.S)
  nt <- length(tseq); t.by <- tseq[2] - tseq[1]
  n.S <- sapply(y.S, length); N <- sum(n.S)
  scoreY <- 0
  for (l in 1:L) {
    xl <- x.S[[l]]; yl <- y.S[[l]]; nl <- n.S[l]
    score <- (xl-matrix(colMeans(xl), nl, nt, byrow=T))%*%evec.est*t.by
    scoreY <- scoreY + nl/N * crossprod(score, yl-mean(yl)) / (nl-1)
  }
  w.est <- scoreY / eval.est
  return(list(w.est=w.est, evec.est=evec.est, eval.est=eval.est))
}

tlFLR_step2 <- function(x, y, tseq, m, lambda, w.est, evec.est) {
  # step 2: de-bias
  # the truncation number m is a scalar
  
  n <- length(y)
  nt <- length(tseq); t.by <- tseq[2] - tseq[1]
  w.est <- w.est[1:m]
  evec.est <- evec.est[, 1:m, drop=F]
  score <- (x-matrix(colMeans(x), n, nt, byrow=T))%*%evec.est*t.by
  yy <- y - mean(y) - score %*% w.est
  if (m==1) {
    b.est <- w.est + sum(score*yy)/sum(score^2)
    slope.est <- evec.est %*% b.est
    lambda <- 0
  } else {
    result.lasso <- glmnet(x=score, y=yy, family='gaussian', lambda=lambda, intercept=FALSE)
    lamlen <- length(lambda)
    b.est <- matrix(0, m, lamlen)
    slope.est <- matrix(0, nt, lamlen)
    for (j in 1:lamlen) {
      b.est[,j] <- w.est + result.lasso$beta[,j]
      slope.est[,j] <- evec.est %*% b.est[,j]
    }
  }
  return(list(slope.est=slope.est, b.est=b.est, evec.est=evec.est, m=m, lambda=lambda))
}

FLR <- function(x, y, x.S=NULL, y.S=NULL, tseq, m){
  
  L <- length(x.S)
  x.all <- vector(mode='list', length = L+1)
  y.all <- vector(mode='list', length = L+1)
  for (l in 1:(L+1)) {
    if (l==1) {
      x.all[[l]] <- x
      y.all[[l]] <- y
    } else {
      x.all[[l]] <- x.S[[l-1]]
      y.all[[l]] <- y.S[[l-1]]
    }
  }
  if (is.null(m)){
    result.FPCA <- estFPCA(x.all, tseq)
  }else{
    result.FPCA <- estFPCA(x.all, tseq, m=max(m))
  }
  evec.est <- result.FPCA$evec.est
  eval.est <- result.FPCA$eval.est
  
  nt <- length(tseq); t.by <- tseq[2] - tseq[1] 
  n.all <- sapply(x.all, nrow); N <- sum(n.all)
  scoreY <- 0
  for (l in 1:(L+1)) {
    xl <- x.all[[l]]; yl <- y.all[[l]]; nl <- n.all[l]
    score.l <- (xl-matrix(colMeans(xl), nl, nt, byrow=T))%*%evec.est*t.by
    scoreY <- scoreY + nl/N * crossprod(score.l, yl-mean(yl)) / (nl-1)
  }
  b.est <- scoreY / eval.est
  slope.est <- evec.est %*% b.est
  return(list(slope.est=slope.est, b.est=b.est, evec.est=evec.est))
}

cv.tlFLR <- function(x, y, w.est, evec.est, tseq, m, lambda, k=5) {
  # k-fold cross-validation for tuning parameters
  # minimize the out-of-sample prediction error using the target dataset
  
  nt <- length(tseq); t.by <- tseq[2] - tseq[1]
  n <- length(y)
  ntest <- floor(n/k); ntrain <- n - ntest
  mlen <- length(m); lamlen <- length(lambda)
  error <- array(0, dim=c(k, mlen, lamlen))
  for (i in 1:k) {
    xtest <- x[((i-1)*ntest+1):(i*ntest),,drop=F]; xtrain <- x[-(((i-1)*ntest+1):(i*ntest)),,drop=F]
    ytest <- y[((i-1)*ntest+1):(i*ntest)]; ytrain <- y[-(((i-1)*ntest+1):(i*ntest))]
    scoretest <- xtest%*%evec.est*t.by
    for (j in 1:mlen) {
      m.tmp <- m[j]
      if (m.tmp==1) {
        result <- tlFLR_step2(xtrain, ytrain, tseq, m.tmp, 0, w.est, evec.est)
        slope.est <- result$slope.est
        b.est <- result$b.est
        a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
        error[i,j,] <- mean((ytest - a.est - scoretest[,1:m.tmp,drop=F]%*%b.est)^2) # out-of-sample prediction error
      }else {
        result <- tlFLR_step2(xtrain, ytrain, tseq, m.tmp, lambda, w.est, evec.est)
        for (l in 1:lamlen) {
          slope.est <- result$slope.est[,l]
          b.est <- result$b.est[,l]
          a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
          error[i,j,l] <- mean((ytest - a.est - scoretest[,1:m.tmp,drop=F]%*%b.est)^2)
        }
      }  
    }
  }
  error_mean <- matrix(0, mlen, lamlen)
  for(i in 1:k){
    error_mean <- error_mean + error[i,,]/k
  }
  idx_lam <- ceiling(which.min(error_mean)/mlen)
  lam.opt <- lambda[idx_lam]
  idx_m <- which.min(error_mean)-(idx_lam-1)*mlen
  m.opt <- m[idx_m]
  if (m.opt==1) {
    lam.opt <- 0
  }
  
  return(list(m.opt=m.opt, lam.opt=lam.opt, error_mean=error_mean, error=error, m=m, lambda=lambda))
}

cv.FLR <- function(x, y, x.S=NULL, y.S=NULL, tseq, m, k=5) {
  # k-fold cross-validation for tuning parameters
  # minimize the out-of-sample prediction error using the target dataset
  
  n <- length(y)
  nt <- length(tseq); t.by <- tseq[2] - tseq[1]
  ntest <- floor(n/k); ntrain <- n - ntest
  mlen <- length(m)
  error <- matrix(0, k, mlen)
  for (i in 1:k) {
    xtest <- x[((i-1)*ntest+1):(i*ntest),,drop=F]; xtrain <- x[-(((i-1)*ntest+1):(i*ntest)),,drop=F]
    ytest <- y[((i-1)*ntest+1):(i*ntest)]; ytrain <- y[-(((i-1)*ntest+1):(i*ntest))]
    result <- FLR(xtrain, ytrain, x.S, y.S, tseq, m=max(m))
    b.est <- result$b.est
    evec.est <- result$evec.est 
    scoretest <- xtest%*%evec.est*t.by
    for(j in 1:mlen){
      m.tmp <- m[j]
      slope.est <- evec.est[,1:m.tmp,drop=F] %*% b.est[1:m.tmp]
      a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
      error[i,j] <- mean((ytest - a.est - scoretest[,1:m.tmp,drop=F]%*%b.est[1:m.tmp])^2)
    }
  }
  error_mean <- colMeans(error)
  idx_min <- which.min(error_mean)
  m.opt <- m[idx_min]
  return(list(m.opt = m.opt, error_mean=error_mean, error=error, m=m))
}

# adaptive algorithm when the informative datasets are unknown

getStat_gfun_error <- function(x, y, x.S, y.S, tseq) {
  # squared L2 norm between the functions of (sample) covariance between X(t) and Y
  t.by <- tseq[2] - tseq[1]
  L <- length(y.S)
  n <- length(y)
  g_est <- (crossprod(x, y-mean(y)) / n)
  error_g <- rep(NA, L)
  for (l in 1:L) {
    xl <- x.S[[l]]; yl <- y.S[[l]]; nl <- length(yl)
    error_g[l] <- sum((crossprod(xl, yl-mean(yl)) / nl - g_est)^2) * t.by
  }
  return(error_g)
}

tlFLR_Agg <- function(x, y, x.S, y.S, tseq, m, lambda, opts, num_splits=5, split_ratio=0.5, k=5, targetInc=TRUE) {
  
  n <- nrow(x); nt <- length(tseq)
  n1 <- ceiling(n*split_ratio)
  n2 <- n - n1 # used for aggregation
  slope.est.qagg <- 0
  slope.est.sagg <- 0
  for(j in 1:num_splits) {
    idx <- sample(1:n, n1, replace=FALSE)
    x1 <- x[idx,]; y1 <- y[idx]
    x2 <- x[-idx,]; y2 <- y[-idx]
    
    # construct candidate sets of transferrable source data
    error <- getStat_gfun_error(x1, y1, x.S, y.S, tseq)
    
    error2 <- unique(sort(error)) # ascending
    n_sets <- length(error2)
    set_list <- list()
    for (l in 1:n_sets) {
      set_list[[l]] <- which(error <= error2[l])
    }
    set_list <- unique(set_list)
    
    # estimates
    slope.est.all <- matrix(0, nt, length(set_list)+1)
    cv_result <- cv.FLR(x1, y1, NULL, NULL, tseq, m, k=k)
    result <- FLR(x1, y1, NULL, NULL, tseq, cv_result$m.opt)
    slope.est.all[,1] <- result$slope.est
    for (l in 1:length(set_list)) {
      idx_set <- set_list[[l]]
      result_step1 <- tlFLR_step1(x1, y1, x.S[idx_set], y.S[idx_set], tseq, max(m), targetInc)
      cv_result <- cv.tlFLR(x1, y1, result_step1$w.est, result_step1$evec.est, tseq, m, lambda, k=k)
      # result_debias <- tlFLR_step2(x1, y1, tseq, cv_result$m.opt, cv_result$lam.opt, result_step1$w.est, result_step1$evec.est)   
      if(cv_result$m.opt==1){
        result_debias <- tlFLR_step2(x1, y1, tseq, cv_result$m.opt, 0, result_step1$w.est, result_step1$evec.est) 
        slope.est.all[,l+1] <- result_debias$slope.est
      } else {
        result_debias <- tlFLR_step2(x1, y1, tseq, cv_result$m.opt, lambda, result_step1$w.est, result_step1$evec.est) 
        idx_lam <- which(lambda == cv_result$lam.opt)
        slope.est.all[,l+1] <- result_debias$slope.est[,idx_lam]
      }
    }
    # Q aggregation
    nparam <- length(opts$omega)
    if (nparam==1) {
      result_qagg <- Qagg(x2, y2, slope.est.all, tseq, opts)
    }else {
      opts_work <- opts
      t.by <- tseq[2] - tseq[1]
      error <- matrix(NA, k, nparam)
      for (idx_k in 1:k) {
        ntest <- floor(n2/k); ntrain <- n2 - ntest
        xtest <- x2[((idx_k-1)*ntest+1):(idx_k*ntest),,drop=F]; xtrain <- x2[-(((idx_k-1)*ntest+1):(idx_k*ntest)),,drop=F]
        ytest <- y2[((idx_k-1)*ntest+1):(idx_k*ntest)]; ytrain <- y2[-(((idx_k-1)*ntest+1):(idx_k*ntest))]
        for(idx_omega in 1:nparam){
          opts_work$omega <- opts$omega[idx_omega]
          slope.est <- Qagg(xtrain, ytrain, slope.est.all, tseq, opts_work)$slope.est
          a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
          error[idx_k, idx_omega] <- mean((ytest - a.est - xtest%*%slope.est*t.by)^2)
        }
      }
      opts_work$omega <- opts$omega[which.min(colMeans(error))]
      result_qagg <- Qagg(x2, y2, slope.est.all, tseq, opts_work)
    }
    slope.est.qagg <- slope.est.qagg + result_qagg$slope.est
    # weight.qagg <- result_qagg$weight
    
    # sparse aggregation
    result_sagg <- sparseAgg(x2, y2, slope.est.all, tseq)
    slope.est.sagg <- slope.est.sagg + result_sagg$slope.est
    # weight.sagg <- result_sagg$weight
  }
  slope.est.qagg <- slope.est.qagg / num_splits
  slope.est.sagg <- slope.est.sagg / num_splits
  
  return(list(slope.est.qagg = slope.est.qagg, slope.est.sagg=slope.est.sagg))
}


FLR_Agg <- function(x, y, x.S, y.S, tseq, m, opts, num_splits=5, split_ratio=0.5, k=5) {
  
  n <- nrow(x); nt <- length(tseq)
  n1 <- ceiling(n*split_ratio)
  n2 <- n - n1 # used for aggregation
  slope.est.qagg <- 0
  slope.est.sagg <- 0
  for(j in 1:num_splits) {
    idx <- sample(1:n, n1, replace=FALSE)
    x1 <- x[idx,]; y1 <- y[idx]
    x2 <- x[-idx,]; y2 <- y[-idx]
    
    # construct candidate sets of transferrable source data
    error <- getStat_gfun_error(x1, y1, x.S, y.S, tseq)
    
    error2 <- unique(sort(error)) # ascending
    n_sets <- length(error2)
    set_list <- list()
    for (l in 1:n_sets) {
      set_list[[l]] <- which(error <= error2[l])
    }
    set_list <- unique(set_list)
    
    # estimates
    slope.est.all <- matrix(0, nt, length(set_list)+1)
    cv_result <- cv.FLR(x1, y1, NULL, NULL, tseq, m, k=k)
    result <- FLR(x1, y1, NULL, NULL, tseq, cv_result$m.opt)
    slope.est.all[,1] <- result$slope.est
    for (l in 1:length(set_list)) {
      idx_set <- set_list[[l]]
      cv_result <- cv.FLR(x1, y1, x.S[idx_set], y.S[idx_set], tseq, m, k=k)
      result <- FLR(x1, y1, x.S[idx_set], y.S[idx_set], tseq, cv_result$m.opt)
      slope.est.all[,l+1] <- result$slope.est
    }
    # Q aggregation
    nparam <- length(opts$omega)
    if (nparam==1) {
      result_qagg <- Qagg(x2, y2, slope.est.all, tseq, opts)
    }else {
      opts_work <- opts
      t.by <- tseq[2] - tseq[1]
      error <- matrix(NA, k, nparam)
      for (idx_k in 1:k) {
        ntest <- floor(n2/k); ntrain <- n2 - ntest
        xtest <- x2[((idx_k-1)*ntest+1):(idx_k*ntest),,drop=F]; xtrain <- x2[-(((idx_k-1)*ntest+1):(idx_k*ntest)),,drop=F]
        ytest <- y2[((idx_k-1)*ntest+1):(idx_k*ntest)]; ytrain <- y2[-(((idx_k-1)*ntest+1):(idx_k*ntest))]
        for(idx_omega in 1:nparam){
          opts_work$omega <- opts$omega[idx_omega]
          slope.est <- Qagg(xtrain, ytrain, slope.est.all, tseq, opts_work)$slope.est
          a.est <- mean(ytrain) - sum(colMeans(xtrain)*slope.est*t.by)
          error[idx_k, idx_omega] <- mean((ytest - a.est - xtest%*%slope.est*t.by)^2)
        }
      }
      opts_work$omega <- opts$omega[which.min(colMeans(error))]
      result_qagg <- Qagg(x2, y2, slope.est.all, tseq, opts_work)
    }
    slope.est.qagg <- slope.est.qagg + result_qagg$slope.est
    # weight.qagg <- result_qagg$weight
    
    # sparse aggregation
    result_sagg <- sparseAgg(x2, y2, slope.est.all, tseq)
    slope.est.sagg <- slope.est.sagg + result_sagg$slope.est
    # weight.sagg <- result_sagg$weight
  }
  slope.est.qagg <- slope.est.qagg / num_splits
  slope.est.sagg <- slope.est.sagg / num_splits
  
  return(list(slope.est.qagg = slope.est.qagg, slope.est.sagg=slope.est.sagg))
}

Qagg <- function(x, y, slope.est.all, tseq, opts) {
  # @param slope.est.all, matrix of all estimators, nt*(L+1)
  # @param opts, list of optimization parameters

  omega <- opts$omega
  niter <- opts$niter
  tol <- opts$tol
  t.by <- tseq[2] - tseq[1]
  L <- ncol(slope.est.all); n <- nrow(x)

  # exponential weight as initial value, temperature parameter: omega
  # smaller omega, sparser weights; larger omega, more uniform weights
  mse.all <- colSums((matrix(y, n, L, byrow = FALSE) - x %*% slope.est.all * t.by)^2)
  weight <- exp(-mse.all/omega) / sum(exp(-mse.all/omega))
  slope.est <- slope.est.all %*% weight
  weight.old <- weight
  iter <- 0; error <- Inf
  while ((iter < niter) & (error > tol)) {
    weight <- exp(-mse.all/omega + 1/(4*omega)*colSums((x %*% slope.est.all * t.by - matrix(x %*% slope.est * t.by, n, L))^2))
    weight <- weight / sum(weight)
    error <- sum((weight-weight.old)^2)
    weight.old <- weight
    # update
    slope.est <- 0.75*slope.est + 0.25 * slope.est.all %*% weight
    iter <- iter + 1
  }
  return(list(slope.est=slope.est, weight=weight))
}

sparseAgg <- function(x, y, slope.est.all, tseq) {
  # without pre-selection step
  # in the section "simulation study" of the original paper
  
  t.by <- tseq[2] - tseq[1]
  n <- nrow(x)
  L <- ncol(slope.est.all)
  mse.all <- colSums((y - x %*% slope.est.all * t.by)^2)/n
  idx.opt <- which.min(mse.all)
  error <- rep(Inf, L)
  for (l in 1:L) {
    if (l==idx.opt) {
      weight <- 1
      error[l] <- sum((y - x %*% slope.est.all[,idx.opt] * t.by)^2) / n
      if (error[l] <= min(error)) {
        slope.est.opt <- slope.est.all[, idx.opt]
        weight.opt <- rep(0, L)
        weight.opt[idx.opt] <- weight
      }
    } else {
      diff <- sum((x %*% slope.est.all[, l] * t.by - x %*% slope.est.all[, idx.opt] * t.by)^2)/n
      weight <- 0.5 * ((mse.all[l] - mse.all[idx.opt])/ diff + 1)
      weight <- min(max(0, weight), 1)
      slope.est <- slope.est.all[,idx.opt]*weight + slope.est.all[,l]*(1-weight)
      error[l] <- sum((y - x %*% slope.est * t.by)^2) / n
      if (error[l] <= min(error)) {
        slope.est.opt <- slope.est
        weight.opt <- rep(0, L)
        weight.opt[idx.opt] <- weight; weight.opt[l] <- 1-weight
      }
    }
  }
  return(list(slope.est=slope.est.opt, weight=weight.opt))
}
