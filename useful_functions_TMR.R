

adj_gen <- function(prob,p){
  adjacency_matrix <- matrix(rbinom(p*p, size = 1, prob = prob), nrow = p, ncol = p)
  return(adjacency_matrix)
}

prec_from_adj <- function(A){
  E <- A
  nb_edges <- sum(E == 1)
  vec_magnitude <- c(-1,1)
  
  E[A == 1] <- runif(nb_edges, min = vec_magnitude[1], max = vec_magnitude[2])
  
  E_bar <- (E + t(E)) / 2
  
  msign <- matrix(1, nrow = nrow(E), ncol = ncol(E))
  msign[upper.tri(msign)] <- sample(c(-1,1), size = sum(upper.tri(msign)),  prob = c(0.5, 0.5), replace = TRUE)
  msign[lower.tri(msign)] <- t(msign)[lower.tri(msign)]
  E_bar <- E_bar * msign
  
  # minimum eigenvalue
  min_eigen <- min(eigen(E_bar, only.values = TRUE)$values)
  if (min_eigen < 0) {
    Omega <- E_bar + (0.1 - min_eigen) * diag(p)
  } else{
    Omega <- E_bar + 0.1 * diag(p)
  }
  
  return(Omega)
  
}

precision = function (g, g.hat) {
  # The proportion of predicted edges that are true
  conf.mat = confusion.matrix(g, g.hat)
  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else if (conf.mat[1, 1] == 0 & conf.mat[1, 2] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[1, 2]))
  }
}

recall = function (g, g.hat) {
  # The proportion of true edges that were identified by the estimated graph
  conf.mat = confusion.matrix(g, g.hat)
  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[2, 1]))
  }
}

GM_gen <-function(n,p,list_hyper, list_init, thresh = 0.02){
  omega <- matrix(0,nrow = p, ncol = p)
  
  lambda <- list_hyper$lambda
  
  a <- list_hyper$a
  b <- list_hyper$b
  
  ar <- list_hyper$ar
  br <-list_hyper$br
  
  v0 <- list_hyper$v0
  v1 <- list_hyper$v1
  
  tau <- rgamma(1,a,b)
  rho <- thresh
  
  var1 <- v1^2/tau
  var0 <- v0^2/tau
  
  for (i in 1:p){
    omega[i,i] <- rexp(1,lambda/2)
    for (j in 1:i-1){
      if (rbinom(1, size =1, prob =  rho) == 1){
        omega[i,j] <- rnorm(1,0,var0)
        omega[j,i] <- omega[i,j]
      } else {
        omega[i,j] <- rnorm(1,0,var1)
        omega[j,i] <- omega[i,j]
      }
    }
  }
  
  min_eigen <- min(eigen(omega, only.values = TRUE)$values)
  if (min_eigen < 0) {
    omega <- omega + (0.1 - min_eigen) * diag(p)
  } else{
    omega <- omega + 0.1 * diag(p)
  }
  
  
  sigma <- solve(omega)
  
  return(list(Sigma = sigma, Omega = omega))
}




GMx_gen <- function(n,p,list_hyper_x, X = list(1,2,3,4)){
  
  if (!is.numeric(n) || !is.numeric(p)) {
    stop("Input arguments must be numeric.")
  }
  
  a <- list_hyper_x$a
  b <- list_hyper_x$br
  
  ar <- list_hyper_x$ar
  br <-list_hyper_x$br
  
  v0 <- list_hyper_x$v0
  v1 <- list_hyper_x$v1
  
  tau <- rgamma(1,a,b)
  
  var1 <- v1^2/tau
  var0 <- v0^2/tau
  
  asig <- list_hyper_x$asig
  bsig <- list_hyper_x$bsig
  
  n0 <- list_hyper_x$n0
  t0 <- list_hyper_x$t0
  
  zeta <- rnorm(1,mean = n0, sd = t0)
  cat('$$$',zeta,'$$$')
  
  sig_neg2 <- rgamma(1, asig, bsig)
  sig <- sig_neg2^-0.5
  
  BETA <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      BETA[i,j] <- rnorm(1,0,sig)
    }
  }
  
  
  RHO <- array(0 , dim = c(p,p,length(X)))
  for (k in 1:length(X)){
    a.ord <- X[[k]]
    for (i in 1:p){
      for (j in 1:i){
        RHO[i,j,k] <- pnorm(zeta + a.ord*BETA[i,j])
        if (i!=j){
          RHO[j,i,k] <- RHO[i,j,k]  
        }
      }
    }
  }
  
  DELTA <- array(0 , dim = c(p,p,length(X)))
  OMEGA <- array(0 , dim = c(p,p,length(X)))
  
  for (k in 1:length(X)){
    for (i in 1:p){
      for (j in 1:i){
        DELTA[i,j,k] <- rbinom(1, size =1, prob = RHO[i,j,k]) 
      }
    }
    for (i in 1:p){
      if (i ==j){
        OMEGA[i,i,k] <- rexp(1,lambda/2)
      }
      else{
        for (j in 1:i){
          if (DELTA[i,j,k] == 1){
            OMEGA[i,j,k] <- rnorm(1,0,var1)
            OMEGA[j,i,k] <- OMEGA[i,j,k]
          } else {
            OMEGA[i,j,k] <- rnorm(1,0,var0)
            OMEGA[j,i,k] <- OMEGA[i,j,k]
          }
        }
      }
    }
    
  }
  
  for (k in 1:length(X)){
    
    
    omega_bar <- (OMEGA[,,k] + t(OMEGA[,,k])) / 2
    
    min_eigen <- min(eigen(omega_bar, only.values = TRUE)$values)
    if (min_eigen < 0) {
      OMEGA[,,k] <- omega_bar + (0.1 - min_eigen) * diag(p)
    } else{
      OMEGA[,,k] <- omega_bar + 0.1 * diag(p)
    }
    
  }
  
  SIGMA = array(0 , dim = c(p,p,length(X)))
  for (k in length(X)){
    SIGMA[,,k] <- solve(OMEGA[,,k])
  }
  
  return(list(Sigma = SIGMA, Omega = OMEGA))
}


perform_ROC_simulation = function(omega.true, n,list_hyper, list_init, N=100, include.glasso=T, include.ssl=T, scale.data = T){
  res=list()
  p=nrow(omega.true)
  sigma.true = solve(omega.true)
  if(include.ssl){
    n.ppi.thresh = 20 # Threshold the posterior inclusion probability to get different FDR (i.e. not thresholding matrix elements)
    ppi.thresh = seq(0,1,n.ppi.thresh)
    res$precisions.ssl = matrix(0,N,n.ppi.thresh)
    res$recalls.ssl = matrix(0,N,n.ppi.thresh)
    res$TPR.ssl = matrix(0,N,n.ppi.thresh)
    res$FPR.ssl = matrix(0,N,n.ppi.thresh)
  }
  if(include.glasso){
    n.lambda = 20 # Varying the penalty parameter to get different FDRs
    res$precisions.glasso = matrix(0,N,n.lambda)
    res$recalls.glasso = matrix(0,N,n.lambda)
    res$TPR.glasso = matrix(0,N,n.lambda)
    res$FPR.glasso = matrix(0,N,n.lambda)
  }
  
  if(include.glasso){
    for(i in 1:N){
      y = mvtnorm::rmvnorm(n, rep(0, p), sigma.true)
      if(scale.data){
        y = scale(y)
      }
      res.glasso = huge(y, method = 'glasso', nlambda=n.lambda, verbose = F)
      res.glasso.omegas = res.glasso$icov # A list of precision matrices
      res$precisions.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$recalls.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$TPR.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      res$FPR.glasso[i,] = unlist(lapply(res.glasso.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
    }
  }
  if(include.ssl){
    for(i in 1:N){
      print(i)
      #res.ssl.omegas <- list()
      #for(x in 1:n.ppi.thresh){
      
      #GM <- GM_gen(n,p,list_hyper, list_init, thresh = x/n.ppi.thresh)
      #adj_mat <- adj_gen(x/(n.ppi.thresh),p)
      adj_mat <- adj_gen(0.08,p)
      new.omega.true <- prec_from_adj(adj_mat)
      new.sigma.true <- solve(new.omega.true)
      #new.omega.true <- GM$Omega
      y = mvtnorm::rmvnorm(n, rep(0, p), new.sigma.true)
      if(scale.data){
        y = scale(y)
      }
      
      res.ssl <- GM(y, list_hyper, list_init) # Must specify the hyper params as well
      res.ssl.omega <- res.ssl$Omega
      #res.ssl.omegas[[x]] <- res.ssl.omega
      
      for(x in 1:n.ppi.thresh){
        res$precisions.ssl[i,x] <- precision(abs(res.ssl.omega)>(x/(2*n.ppi.thresh)), abs(new.omega.true)>(x/(2*n.ppi.thresh)))
        res$recalls.ssl[i,x] <- recall(abs(res.ssl.omega)>(x/(2*n.ppi.thresh)), abs(new.omega.true)>(x/(2*n.ppi.thresh)))
        res$TPR.ssl[i,x] <- TPR(abs(res.ssl.omega)>(x/(2*n.ppi.thresh)), abs(new.omega.true)>(x/(2*n.ppi.thresh)))
        res$FPR.ssl[i,x] <- FPR(abs(res.ssl.omega)>(x/(2*n.ppi.thresh)), abs(new.omega.true)>(x/(2*n.ppi.thresh)))
        
        #res.ssl.omegas <- lapply(ppi.thresh, FUN = function(s) omegas.list$PPI < s) # A list of prec.matrices. Also probably not how you get the PPI, fix this
        #res$precisions.ssl[i,x] = unlist(lapply(res.ssl.omegas, FUN = function(s) precision(abs(s)>1e-5, abs(omega.true)>1e-5)))
        #res$recalls.ssl[i,x] = unlist(lapply(res.ssl.omegas, FUN = function(s) recall(abs(s)>1e-5, abs(omega.true)>1e-5)))
        #res$TPR.ssl[i,x] = unlist(lapply(res.ssl.omegas, FUN = function(s) TPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
        #res$FPR.ssl[i,x] = unlist(lapply(res.ssl.omegas, FUN = function(s) FPR(abs(s)>1e-5, abs(omega.true)>1e-5)))
      }
    }
  }
  
  
  # Average over all N simulations (i.e., find the average performance for each lambda/PPI)
  # Glasso
  if(include.glasso){
    res$mean.precisions.glasso = colMeans(res$precisions.glasso)
    res$mean.recalls.glasso = colMeans(res$recalls.glasso)
    res$mean.TPR.glasso = colMeans(res$TPR.glasso)
    res$mean.FPR.glasso = colMeans(res$FPR.glasso) 
  }
  # Spike-and-slab
  if(include.ssl){
    res$mean.precisions.ssl = colMeans(res$precisions.ssl)
    res$mean.recalls.ssl = colMeans(res$recalls.ssl)
    res$mean.TPR.ssl = colMeans(res$TPR.ssl)
    res$mean.FPR.ssl = colMeans(res$FPR.ssl) 
  }
  return(res)
}