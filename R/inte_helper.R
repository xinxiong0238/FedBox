#==============================================================================
# Notification:
# 1. Each site must have the same number of d_phi and d_psi
# 2.Variable order: psi, phi, X, X*phi
#==============================================================================

pre_inte_fromLocal <- function(local_lst){
  H_lst <- vector('list', 1)
  g_lst <- vector('list', 1)
  theta_lst <- vector('list', 1)
  K <- length(local_lst)
  indx_country <- c()
  n_lst <- c()
  for (k in 1:K) {
    H_lst[[k]] <- local_lst[[k]]$Hessian
    g_lst[[k]] <- local_lst[[k]]$grad
    theta_lst[[k]] <- c(local_lst[[k]]$eta, local_lst[[k]]$xi, local_lst[[k]]$beta)
    indx_country <- c(indx_country, k)
    n_lst <- c(n_lst, nrow(local_lst[[k]]$dat.obs))
  }
  
  
  d_psi <- length(local_lst[[1]]$eta)
  d_phi <- length(local_lst[[1]]$xi)
  p <- length(g_lst[[1]])
  penalty.factor <- c(rep(0, d_psi + d_phi), rep(1, p - d_psi - d_phi))
  names(penalty.factor) <- names(sum_stat$grad)
  
  # Specify on which variables (the main effects) we impose they are the same across the sites.
  feature_name <- names(g_lst[[1]])[!str_detect(names(g_lst[[1]]), "phi|psi")]
  share_var <- c(rep(0, d_psi + d_phi), rep(1, length(feature_name)), 
                 rep(0, p - d_psi - d_phi - length(feature_name)))
  names(share_var) <- names(sum_stat$grad)
  
  return(list(`g_lst` = g_lst, `H_lst` = H_lst, `theta_lst` = theta_lst,
              `n_lst` = n_lst, `share_var` = share_var, `penalty` = penalty.factor))
}

pre_inte_XY <- function(H_lst, g_lst, theta_lst, n_lst, 
                        share_var, penalty.factor){
  M <- length(g_lst)
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  options(warn = -1)
  n <- sum(n_lst)
  p <- length(H_lst[[1]][1, ])
  
  for (m in 1:M){
    H <- H_lst[[m]]
    g <- - g_lst[[m]] + H_lst[[m]] %*% theta_lst[[m]]
    n_m <- n_lst[m]
    
    mat_all <- cbind(H, g)
    mat_all <- rbind(mat_all, t(c(g, max(H) + 10)))
    mat_all <- n_m * mat_all
    
    svd_result <- svd(mat_all)
    s_value <- svd_result$d
    s_mat <- diag(sqrt(s_value))[,1:(1 + min(p, n_lst[m]))]
    data_all <- svd_result$u %*% s_mat
    X <- t(data_all[-length(data_all[ ,1]),])
    Y <- data_all[length(data_all[ ,1]),]
    X_lst[[m]] <- X
    Y_lst[[m]] <- Y
  }
  
  X_all_left <- c()
  X_all_right <- X_lst[[1]][,which(share_var == 0)]
  Y_all <- c()
  penalty.factor.all <- penalty.factor[which(share_var == 1)]
  for (m in 1:M){
    Y_all <- c(Y_all, Y_lst[[m]])
    if (m >= 2){
      X_all_right <- bdiag(X_all_right, X_lst[[m]][,which(share_var == 0)])
    }
    X_all_left <- rbind(X_all_left, X_lst[[m]][,which(share_var == 1)])
    penalty.factor.all <- c(penalty.factor.all, penalty.factor[which(share_var == 0)])
  }
  X_all <- as.matrix(cbind(X_all_left, X_all_right))
  
  return(list(`X_all` = X_all, `Y_all` = Y_all))
}



integrative_fit_for_cv <- function(X_all, Y_all, M, share_var, penalty.factor.lambda, 
                                   lambda1 = 0.1, lambda_init_beta = 0, gamma = 0.5,
                                   max.iter = 10, tol = 1e-5){
  
  # Initial estimator:
  
  fit.init <- glmnet(X_all, Y_all, intercept = F, standardize = F, lambda = 0)
  beta.init <- fit.init$beta
  beta.fit <- beta.init
  
  error_delta <- Inf 
  iter <- 0
  tau <- gamma^(gamma / (1 - gamma)) * (1 - gamma)
  
  while (error_delta > tol & iter < max.iter) {
    iter <- iter + 1
    beta.prev <- beta.fit
    
    theta_mat <- c()
    for (m in 1:M) {
      theta_m <- rep(0, length(share_var))
      theta_m[which(share_var == 1)] <- beta.fit[1:length(which(share_var == 1))]
      indx_m_start <- length(which(share_var == 1)) + 1 + (m - 1) * length(which(share_var == 0))
      indx_m_end <- length(which(share_var == 1)) + m * length(which(share_var == 0)) 
      theta_m[which(share_var == 0)] <- beta.fit[indx_m_start:indx_m_end]
      theta_mat <- cbind(theta_mat, theta_m)
    }
    
    theta_l1_norm <- rowMeans(abs(theta_mat))[which(share_var == 0)]
    
    # Update beta
    
    zeta_vec <- rep(((1 - gamma) / (tau * gamma))^(gamma) * theta_l1_norm^(gamma), M)
    share_factor <- 1 / abs(beta.init[1:length(which(share_var == 1))])
    share_factor <- share_factor / mean(share_factor)
    pen_factor_vec <- c(share_factor, zeta_vec^(1 - 1 / gamma))  ##c(adap lasso weight, groupbidge)
    
    
    ## psi and phi have no penalty, X has extra penalty (lambda2)
    ## *penalty.factor - > c(adap lasso penalty, group bridge)
    pen_vec <- pen_factor_vec * penalty.factor.lambda 
    indx_use <- which(pen_vec != Inf)
    pen_use <- pen_vec[indx_use]
    fit.m <- glmnet(X_all[,indx_use], Y_all, intercept = F, standardize = F, 
                    lambda = lambda1 / nrow(X_all) * mean(pen_use),
                    alpha = 1, penalty.factor = pen_use)
    beta.fit <- rep(0, ncol(X_all))
    beta.fit[indx_use] <- fit.m$beta
    
    # Error evaluation
    error_delta <- (norm(beta.fit - beta.prev, type = '2'))^2
    #print(error_delta)
  }
  
  theta_mat <- c()
  for (m in 1:M) {
    theta_m <- rep(0, length(share_var))
    names(theta_m) = names(share_var)
    theta_m[which(share_var == 1)] <- beta.fit[1:length(which(share_var == 1))]
    indx_m_start <- length(which(share_var == 1)) + 1 + (m - 1) * length(which(share_var == 0))
    indx_m_end <- length(which(share_var == 1)) + m * length(which(share_var == 0)) 
    
    theta_m[which(share_var == 0)] <- beta.fit[indx_m_start:indx_m_end]
    theta_mat <- cbind(theta_mat, theta_m)
  }
  
  return(theta_mat)
}


derive_local = function(time,event,x,calendar_time, num_bh_knot=2, d_phi=3, d_psi,
                        theta){
  
  phi = data.frame(ns(calendar_time,df=d_phi))
  colnames(phi) = paste0('phi',1:ncol(phi))
  
  # make it shifted to be orthogonal to 1
  for (j in 1:d_phi) {
    phi[,j] <- (phi[,j] - mean(phi[,j])) / sd(phi[,j])
  }
  
  # interaction with covariates
  X_ob0 = x
  X_ob = X_ob0
  for (j in 1:d_phi) {
    inter_basis = phi[,j]*X_ob0
    colnames(inter_basis) = paste0(colnames(phi)[j],".",colnames(X_ob0))
    X_ob = X_ob %>% mutate(inter_basis)
  }
  # interaction with missing indicator
  dat.obs = mutate(cbind(X_ob, phi),y=time,failed=event,calendar_time=calendar_time)
  local_est = get_local_est(dat.obs,num_bh_knot,d_phi)
  ### use estimations to compute summary statistics
  sum_stat = summary_stat(Delta = dat.obs$failed,
                          time = dat.obs$y,
                          X = local_est$X,
                          beta = local_est$beta, eta = local_est$eta,
                          xi = local_est$xi, phi = local_est$phi,
                          psi = local_est$psi, cbhhat.obj = local_est$cbhhat.obj)
  return(c(local_est,sum_stat))
}





Cal_GIC <- function(H_lst, g_lst, theta_lst, n_lst, share_var, theta_mat, penalty.factor, 
                    lambda1 = 0.1, gamma = 0.5, type = 'BIC'){
  M <- length(g_lst)
  options(warn = -1)
  n <- sum(n_lst)
  p <- length(H_lst[[1]][1, ])
  tau <- gamma^(gamma / (1 - gamma)) * (1 - gamma)
  zeta_vec <- ((1 - gamma) / (tau * gamma))^(gamma) * (rowSums(abs(theta_mat)))^(gamma)
  Design_mat_all <- c()
  W_mat_all <- c()
  
  #print(rowSums(abs(theta_mat)))
  
  # Calculate df:
  
  indx_diff <- which(share_var == 0 & penalty.factor > 0)
  Design_mat_all <- n_lst[1] * H_lst[[1]][indx_diff,indx_diff]
  W_mat_all <- diag(lambda1 * (zeta_vec[indx_diff]^(1 - 1 / gamma) + 1e-20) / 
                      (abs(theta_mat[indx_diff,1]) + 1e-40))
  
  for (m in 2:M){
    Design_mat_all <- bdiag(Design_mat_all, 
                            n_lst[m] * H_lst[[m]][indx_diff,indx_diff])
    W_mat_all <- bdiag(W_mat_all, 
                       diag(lambda1 * (zeta_vec[indx_diff]^(1 - 1 / gamma) + 1e-20) /
                              (abs(theta_mat[indx_diff,m]) + 1e-40)))
  }
  
  df <- sum(diag(solve(W_mat_all + Design_mat_all) %*% Design_mat_all)) 
  #print(df)
  df <- df + length(which(zeta_vec > 0 & share_var == 1))
  #print(df)
  GIC <- 0
  for (m in 1:M) {
    H <- H_lst[[m]]
    g <- - g_lst[[m]] + H_lst[[m]] %*% theta_lst[[m]]
    n_m <- n_lst[m]
    GIC <- GIC + n_m * t(theta_mat[,m]) %*% H %*% theta_mat[,m] - 2* n_m * t(theta_mat[,m]) %*% g 
  }
  
  if (type == 'BIC'){
    GIC <- GIC / sum(n_lst) + df * log(sum(n_lst)) / sum(n_lst)
  }
  if (type == 'AIC'){
    GIC <- GIC / sum(n_lst) + df * 2 / sum(n_lst)
  }
  if (type == 'RIC'){
    GIC <- GIC / sum(n_lst) + df * log(p) / sum(n_lst)
  }
  if (type == 'mBIC'){
    GIC <- GIC / sum(n_lst) + df * log(max(log(p), exp(1))) * log(sum(n_lst)) / sum(n_lst)
  }
  return(GIC)
}
