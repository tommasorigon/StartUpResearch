Gibbs_model1 <- function(B,W,prior,R,burn_in,thinning,verbose=TRUE){
  
  verbose_step = ceiling(R/50)
  
  n_t <- ncol(B)
  n_l <- nrow(B)
  n   <- n_l * n_t # Total number of rows and colums
  
  # Not missing values
  N_l <- apply(B,2,function(x) sum(!is.na(x)))
  N   <- sum(N_l)  # Total number of observations
  
  w_plus  <- rowSums(W)
  
  # Hyperparameters
  sigma2_0 <- prior$sigma2_0
  psi      <- prior$psi
  a_sigma2      <- prior$a_sigma2
  b_sigma2      <- prior$b_sigma2
  a_gamma2      <- prior$a_gamma2
  b_gamma2      <- prior$b_gamma2
  a_tau2        <- prior$a_tau2
  b_tau2        <- prior$b_tau2
  
  
  R_corr       <- exp(-psi*as.matrix(dist(1:n_t)))
  P_corr       <- solve(R_corr)
  
  # Raise an error if the matrix inversion seems to be unstable
  if(abs(sum((P_corr%*%R_corr)^2) - n_t) > 1e-10) stop("The covariance matrix of the Gaussian Process is numerically ill-conditioned")
  
  # Initialization
  tau2 <-  sigma2 <- gamma2 <- 1
  phi  <-  numeric(n_l)
  Z    <-  numeric(n_t)
  
  # Output
  beta_0_out <- numeric(R)
  Z_out      <- matrix(0,R,n_t)
  phi_out    <- matrix(0,R,n_l)
  gamma2_out <- numeric(R)
  tau2_out   <- numeric(R)
  sigma2_out <- numeric(R)
  
  # Gibbs sampling
  for(r in 1:(R*thinning + burn_in)){
    
    # First step (Intercept)
    mu_tilde     <- sigma2_0 * sum( t(B - phi) - Z,na.rm = TRUE) / (gamma2 + N*sigma2_0)
    sigma2_tilde <- gamma2*sigma2_0/(gamma2 + N*sigma2_0)
    beta_0       <- rnorm(1, mu_tilde,sqrt(sigma2_tilde))
    
    # Second step  (CAR model)
    for(l in 1:n_l) {
      mu_tilde     <- (tau2/w_plus[l]*sum( B[l,] - beta_0 - Z,na.rm = TRUE) + gamma2*sum(W[l,-l]*phi[-l]/w_plus[l]))/ (gamma2+tau2*n_t/w_plus[l])
      sigma2_tilde <-  (gamma2*tau2/w_plus[l])/(gamma2 + n_t*tau2/w_plus[l])
      phi[l]       <- rnorm(1,mu_tilde,sqrt(sigma2_tilde))
    }
    phi <- phi - mean(phi)
    
    # Third step (Gaussian Process)
    Bl_sum       <- colSums(B - beta_0 - phi,na.rm=TRUE)
    eig         <- eigen(diag(N_l/gamma2) + P_corr/sigma2, symmetric = TRUE)
    Sigma_tilde <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_tilde    <- c(Sigma_tilde %*% (Bl_sum/gamma2))
    Z           <- mu_tilde + c(matrix(rnorm(1 * n_t), 1, n_t) %*% (t(eig$vectors)/sqrt(eig$values)))
    
    # Variances
    tau2   <- 1/rgamma(1,a_tau2 + n_l/2, b_tau2 + sum(phi^2)/2)
    gamma2 <- 1/rgamma(1,a_gamma2 + N/2, b_gamma2 + sum((t(B - beta_0 - phi) - Z)^2,na.rm = TRUE)/2)
    
    mahalanob <- c(t(Z)%*%P_corr%*%Z)
    sigma2    <- 1/rgamma(1,a_sigma2 + n_t/2,b_sigma2 + mahalanob/2)
    
    # Output
    if (r > burn_in & ((r - burn_in) %% thinning == 0)) {
      rr             <- floor((r - burn_in)/thinning)
      beta_0_out[rr] <- beta_0
      Z_out[rr,]     <- Z
      phi_out[rr,]   <- phi
      gamma2_out[rr] <- gamma2
      tau2_out[rr]   <- tau2
      sigma2_out[rr] <- sigma2
    }
    if (verbose) {
      if (r%%(verbose_step*thinning) == 0) 
        cat(paste("Sampling iteration: ", r, " out of ",R*thinning + burn_in, "\n",
                  sep = ""))
    }
  }
  return(list(beta_0 = beta_0_out,
              Z=Z_out,
              phi=phi_out,
              gamma2= gamma2_out,
              tau2=tau2_out,
              sigma2=sigma2_out
  ))
}
