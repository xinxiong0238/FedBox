#================================================================================
#===========   Generate Multiple Local data sets  ===============================
#================================================================================

generate_data_Multilocals <- function(seed = 123, n_vec, T_vec,
                                      beta_nonzero_num, beta_zero_num){
  set.seed(seed)
  K = length(n_vec)
  dist_family <- 'gompertz'
  lambda <- 0.05
  n_beta <- beta_nonzero_num + beta_zero_num
  Xcov <- 1 * ar1_cor(n_beta, 0.2)
  
  betas_name <- c(paste0("X", 1:beta_nonzero_num),
             paste0("M", 1:beta_nonzero_num),
             paste0("gw.X", 1:beta_nonzero_num),
             paste0("gw.M", 1:beta_nonzero_num))
  betas <- data.frame("coeff_name" = betas_name, 
                      "coeff" = c(runif(round(beta_nonzero_num/2), min = -2, max = -0.1),
                                  runif(beta_nonzero_num - round(beta_nonzero_num/2),
                                                                 min = 0.1, max = 2),
                                  runif(beta_nonzero_num, min = -0.3, max = 0.3),
                                  runif(beta_nonzero_num, min = -2, max = 2),
                                  runif(beta_nonzero_num, min = -0.4, max = 0.4)
                                  )
                      )
  rownames(betas) <- betas$coeff_name; betas$coeff_name = NULL
  betas = t(betas)
  sample_data_list <- c()
  for(k in 1:K){
    T0 <- T_vec[k]
    n_patient <- n_vec[k]
    beta_noise <- rnorm(ncol(betas), 0, 0.1)
    beta_noise[1:beta_nonzero_num] = 0
    rho <- runif(1,0.3,0.6)
    local_betas <- beta_noise + betas
    missing_rate <- sample(c(0.1,0.2,0.3), n_beta, replace = TRUE, 
                           prob = c(0.6, 0.3, 0.1))
    
    
    ### simulate data for each site given parameters
    dat_obj <- generate_data(n_patient, n_beta, local_betas,
                             T0, missing_rate, lambda, rho,
                             Xcov, dist_family)
    sample_data <- cbind(dat_obj$dat.obs[, c('id', 'y', 'failed',
                                   paste0('X', 1:n_beta), 'calendar_time')],
                         `lp` = dat_obj$lp.risk)
    sample_data_list <- c(sample_data_list, list(sample_data))
  }
  return(sample_data_list)
}


generate_data = function(n_patient, n_beta, betas, T0,
                         missing_rate, lambda, rho, Xcov,
                         dist_family='weilbull'){
  ### the main function used to generate survival data
  ### Input are parameters 
  ### Output: dat.true: real but unobserved covariates used to generate survival time; dat.obs: observable covariates
  X0 <- mvrnorm(n_patient, mu = rep(0, n_beta), Sigma = Xcov) # mvn with cov=ar(1)
  colnames(X0) <- paste0('X', 1:ncol(X0))
  miss_ind <- matrix(0, nrow = nrow(X0), ncol = ncol(X0)) # missing indicator
  colnames(miss_ind) <- paste0('M', 1:ncol(X0))
  X <- cbind(X0, miss_ind) # combine startum, mvn & missing indicator
  for (j in 1:n_beta){
    miss_indx <- X0[,j] < 1
    X[miss_indx, paste0('M',j)] <- rbinom(sum(miss_indx),1,missing_rate[j])
    X[X[,n_beta+j] == 1, paste0('X',j)] <- 0
  }
  
  calendar_time <- sample(1:T0, n_patient, replace = T) # sample calendar time
  
  gw_value <- gw(calendar_time)
  interact <- gw_value * X # add interaction term
  colnames(interact) <- paste0('gw.', colnames(X))
  
  # interaction with missing indicator
  X_inter <- cbind(X, interact)
  X_true <- cbind(X_inter, data.frame(gw = gw_value)) # combine coviate, missing indicator, interaction and calendar time
  beta_all <- matrix(0, ncol=ncol(X_true))
  colnames(beta_all) <- colnames(X_true)
  for (p in colnames(betas)){
    beta_all[1,p] <- betas[1, p]
  }
  
  # use the weibull distribution to get survival time and status
  dat <- simul(N = n_patient, lambda = lambda, rho = rho, 
               beta = as.vector(beta_all), X = X_true, 
               dist_family = dist_family)
  dat.true <- cbind(dat, calendar_time)
  lp <- as.matrix(X_true) %*% as.matrix(as.vector(beta_all))
  
  # generate observed data
  dat.obs <- cbind(subset(dat, select = c(id, y, failed)), X, calendar_time)
  
  return(list(dat.true = dat.true,
              dat.obs = dat.obs,
              lp.risk = lp))
}

simul <- function(N, X, betas, lambda, rho, beta, dist_family='weibull'){
  ### get Weibull/Gompertz distribution with cumulative baseline hazard: H(t)=lambda*t^rho/H(t)=lambda/rho*(exp(rho*t)-1)
  ### rho changes across strata, for simplicity only assume 2 strata here
  ### formula to simulate time with Weibull, T using uniform distribution, v: (-log(v)/(lambda*exp(beta*X)))^(1/rho)
  ### to simulate time with Gompertz: ((1/rho)*log(1-(rho*log(v))/(lambda*exp(beta*X))))
  Tlat <- rep(0,N)
  cbh <- rep(0,N)
  v <- runif(n=N)
  # Weibull latent time
  if(dist_family=='weibull'){
    Tlat = (-log(v)/(lambda*exp(as.matrix(X)%*%as.matrix(beta))))^(1/rho)
  }
  if(dist_family=='gompertz'){
    Tlat = (1/rho)*
      log(1-(rho*log(v))/
            (lambda*exp(as.matrix(X)%*%as.matrix(beta))))
  }
  # censoring times
  C = runif(n=N,quantile(Tlat,0.5),quantile(Tlat,0.9))
  
  # follow-up times and event indicators
  time <- pmin(Tlat,C)
  status <- as.numeric(Tlat <= C)
  
  # get theoretical value of cumulative baseline hazard and baseline hazard at each time
  return(cbind(data.frame(id=1:N,y=time,failed=status),X))
}

ar1_cor <- function(n, rho) {
  ### construct autoregressive correlation structure
  exponent <- abs(matrix(1:n-1, nrow=n, ncol=n, byrow=TRUE) - (1:n-1))
  rho^exponent
}
gw = function(w){
  ### interacting part of shape function
  cos(w*pi/300)
}


