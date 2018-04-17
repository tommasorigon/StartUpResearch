# The package "klin" perform fast linear algebra involving kronecker products
library(klin)

# --------------------------------------------------------
# - This function allows to obtain the MLE of the PPCA model (Tipping and Bishop, 1999)
# - By default, it returns a rotation for W_ML such that it is lower triangular with positive diagonal (Cholesky decomposition).
# --------------------------------------------------------

# dataset are the raw data
# K is the number of component to be extracted
# If cor=TRUE, the correlation matrix is used instead

PPCA <- function(dataset, K, cor=FALSE){
  
  # Columns of the dataset
  p <- ncol(dataset)
  if(K >= p) stop("K must be strictly lower than the number of variables")
  
  S <- cov(dataset)
  if(cor){S <- cov2cor(S)}
  
  # Eigen-decomposition of the covariance (correlation) matrix
  eigS   <- eigen(S)
  lambda <- pmax(eigS$values,1e-16) # Numerical accuracy
  U      <- eigS$vectors[,1:K]
  
  # Maximum likelihood estimate as in (7) of Tipping and Bishop, (1999)
  sigma2 <- mean(lambda[(K+1):p])                # MLE of the variance
  
  # If K=1, the solution is simply given by. No rotation is needed.
  if(K==1) {
    A = U*sqrt(lambda[1:K] - sigma2)
    if(A[1] < 0){ A <- -1*A}
    return(list(A=t(t(A)),tau=log(sigma2)))
  }
  
  # Otherwise, if K > 1 we get
  W_ML   <- U%*%diag(sqrt(lambda[1:K] - sigma2)) # Using the notation of Tipping and Bishop, we set implicitely  R = Identity
  
  # Rotation to obtain the cholesky decompositionis performed via QR decomposition. 
  # See Property 3.3, Quarteroni et al. (2000)
  A   <- t(qr.R(qr(t(W_ML))))
  A[,diag(A) < 0] <- -1*A[,diag(A) < 0] # The decomposition is adjusted so that the diagonal values are positive
  
  # Notice that Sigma = AA^T = W_ML W_ML^T.
  # The variance in reported in the log-scale
  return(list(A=A,tau=log(sigma2)))
}

# -------------------------------------------------------
# - This function allows to compute standard-errors of the PPCA model
# - It is based on the nonparametric Boostrap (Efron, 1979)
# -------------------------------------------------------

bootPPCA <- function(dataset,K,R=1000,cor=FALSE){
  
  n          <- nrow(dataset) # Number of observation
  indexA     <- lower.tri(matrix(0,ncol(dataset),K),diag=TRUE) 
  param_boot <- matrix(0,R,sum(indexA)+1) # Output, in vectorized form. The logarithm of the varianze is the last element
  
  for(r in 1:R){
    index          <- sample(n,n,replace=TRUE)
    out            <- PPCA(dataset[index,],K)
    param_boot[r,] <- c(out$A[indexA],out$tau)
  }
  
  # The output consist in a matrix containing the bootstrap replicates. 
  param_boot
}

# ------------------------------------------------------------
# - This function compute the log-likelihood of our model, up to an additive constant
# ------------------------------------------------------------

loglikelihood <- function(A,tau,BeigenT,eigenT){
  
  # Eigen decompositions of the matrix Sigma_A
  Sigma_A        <- A%*%t(A)
  eigenA         <- eigen(Sigma_A)
  eigenA$values  <- pmax(eigenA$values,1e-16) # pmax() is necessary to ensure numerical stability
  
  outer_prod <- t(outer(eigenT$values,eigenA$values))
  logdet     <- sum(log(outer_prod + exp(tau)))
  
  x         <- c(crossprod(eigenA$vectors,BeigenT))
  lambdas   <- c(outer_prod) + exp(tau)
  mahalanob <- c(crossprod(x/lambdas,x))
  
  # Log-likelihood
  loglikelihood      <- -0.5*logdet - 0.5*mahalanob
  
  return(loglikelihood)
}

# --------------------------------------------------
# - This function compute the logarithm of the prior distribution of our model, up to an additive constant
# --------------------------------------------------

logprior <- function(A,tau,gamma2,a_sigma,b_sigma,df_diagA){
  
  # Log-prior distribution for the elements of A
  logpriorA   <- sum(dnorm(A[lower.tri(A)],0,sqrt(gamma2),log=TRUE)) + sum(dgamma(diag(A)^2,shape=df_diagA/2,scale=2*gamma2,log=TRUE))
  
  # Log-prior for tau
  logpriortau <- -a_sigma*tau - b_sigma/exp(tau)
  
  # Output
  return(logpriorA + logpriortau)
}

# --------------------------------------------------
# - Main function, which compute the MCMC chain. 
# - Please be patient: this function can be slow, depending on the dimension of the dataset and the power of your computer ;)
# --------------------------------------------------

Metropolis_MCMC <- function(B, time_grid, prior, Iter, burn_in, thinning, tuning, start = NULL, verbose=TRUE){
  
  verbose_step = ceiling(Iter/15) # This shows 15 times the current state of the chain, if verbose=TRUE
  
  # Number of rows and colums
  n_t <- ncol(B)
  n_l <- nrow(B)
  n   <- n_l * n_t # Total number of observations in the dataset
  
  # Hyperparameters
  gamma2      <- prior$gamma2   # Variance of components in A
  psi         <- prior$psi      # Smoothing parameter
  a_sigma     <- prior$a_sigma # Parameters of
  b_sigma     <- prior$b_sigma
  K           <- prior$K
  
  # Correlation matrix for the time
  Sigma_T  <- exp(-psi*as.matrix(dist(time_grid)))
  eigenT   <- eigen(Sigma_T)
  eigenT$values <- pmax(eigenT$values,1e-16)
  
  BeigenT <- B%*%eigenT$vectors

  # If "initialization "start" is not given, use PPCA to initialize the algorithm as default
  if(is.null(start)){start <- PPCA(t(B),K)}
  
  # Initialization and other necessary quantities
  acceptance_ratio      <- 0 # Number of accepted iterations             
  cholSigmaMetropolis   <- tuning$cholSigmaMetropolis # Tuning matrix (Already in the Cholesky form)
  A      <- start$A; A_star <- A      # Initial value for A and A_star
  tau    <- start$tau                 # Initial value for tau
  indexA <- lower.tri(A,diag=TRUE) # Index for the lower triangular part of A
  df_diagA     <- K-1:K+1 # Degrees of freedom for the prior distribution of the diagonal elements of A
  
  # Likelihood and log-posterior
  loglik_old <- loglikelihood(A,tau,BeigenT,eigenT) 
  lp_old     <- loglik_old + logprior(A,tau,gamma2,a_sigma,b_sigma,df_diagA) # Initial value for the log-posterior
  
  # Total number of paramters
  p       <- sum(indexA) + 1 
  
  # Initialization of the output
  sigma2_out <- numeric(Iter)
  A_out      <- array(0,c(Iter,n_l,K))
  logpost    <- numeric(Iter)
  loglik     <- numeric(Iter)
  
  # Metropolis-Hastings for cycle
  for(r in 1:(Iter*thinning + burn_in)){
      
    
      # Sample from the proposal
      sample_proposal <- c(c(A[indexA],tau) + rnorm(p)%*%cholSigmaMetropolis)
      A_star[indexA] <- sample_proposal[-p]
      tau_star       <- sample_proposal[p]
      
      # Computing the new log-posterior
      loglik_new <- loglikelihood(A_star,tau_star,BeigenT,eigenT)
      lp_new     <- loglik_new + logprior(A_star,tau_star,gamma2,a_sigma,b_sigma,df_diagA) # New log-posterior
      alpha      <- exp(lp_new - lp_old) # Acceptance probability
      if(alpha > runif(1)) {A <- A_star; tau <- tau_star; lp_old <- lp_new; loglik_old <- loglik_new; acceptance_ratio <- acceptance_ratio + 1}
      
      # Storing the results
      if (r > burn_in & ((r - burn_in) %% thinning == 0)) {
        rr             <- floor((r - burn_in)/thinning)
        logpost[rr]    <- lp_old
        loglik[rr]     <- loglik_old
        A_out[rr,,]    <- A
        sigma2_out[rr] <- exp(tau)
      }
      # Printing the results
      if (verbose) {
        if (r%%(verbose_step*thinning) == 0) 
          cat(paste("Sampling iteration: ", r, " out of ",Iter*thinning + burn_in, "\n",
                    sep = ""))
      }
  }
  
  # Final Output
  out <- list(logpost = logpost,
              loglik_hat = loglikelihood(apply(A_out,c(2,3),mean),mean(log(sigma2_out)),BeigenT,eigenT),
              loglik  = loglik,
       A = A_out,
       sigma2 = sigma2_out,
       acceptance_ratio = acceptance_ratio/(Iter*thinning + burn_in),
       time_grid = time_grid,
       eigenT = eigenT,
       p = p,
       psi=psi)
  class(out) <- "model3" #Code-name for our model
  return(out)
}

# --------------------------------------------------------
# - Function for computing the Information Criteria DIC
# --------------------------------------------------------

IC <- function(model){
  p      <- model$p
  p_DIC  <- 2*(model$loglik_hat - mean(model$loglik))
  DIC    <- -2*model$loglik_hat + 2*p_DIC
  
  return(cbind(DIC,p,p_DIC))
}

# --------------------------------------------------------------
# - This functions allows to predict the multivariate process, given the MAP
# - If now new time grid is supplied, the function return the prediction over the same set of data
# --------------------------------------------------------------

predict.model3 <- function(object,data,new_grid=NULL){
  
  # The predictions are always based on the MAP
  id_MAP <- which.max(object$logpost)
  A      <- object$A[id_MAP,,]
  sigma2 <- object$sigma2[id_MAP]
  eigenT <- object$eigenT
  B      <- data
  psi    <- object$psi
  
  Sigma_A       <- A%*%t(A)
  eigenA        <- eigen(Sigma_A)
  eigenA$values <- pmax(eigenA$values,1e-16)
  
  # If no new time grid is available
  if(is.null(new_grid)){
    eigenT0 <- eigenT
    n_obs   <- n_new <- length(object$time_grid)
  } else {
    n_obs    <- length(object$time_grid)
    n_new    <- length(new_grid)
    Sigma_T0 <- exp(-psi*as.matrix(dist(c(new_grid, object$time_grid)))[1:n_new, -(1:n_new)])
  }

  x             <- c(t(eigenA$vectors)%*%B%*%eigenT$vectors)
  lambda        <- c(t(outer(eigenT$values,eigenA$values)))

  C1vecB        <- c(eigenA$vectors%*%tcrossprod(matrix(x/(lambda + sigma2), ncol=n_obs),eigenT$vectors))
  pred          <- klin.eval(list(Sigma_T0,Sigma_A),C1vecB) #A fast version of (Sigma_T0 %x% Sigma_A)%*%C1vecB
  return(matrix(pred,ncol=n_new))
}

