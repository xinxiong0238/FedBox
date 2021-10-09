integrative_fit <-  function(inte_data, 
                                lambda1_lst = 0.1, lambda2_lst = lambda1_lst,
                                lambda_init_beta = 0, gamma = 0.5,
                                max.iter = 10, tol = 1e-5, type = 'BIC'){
  ## Take local data
  M <- length(inte_data$n_lst)
  g_lst <- inte_data$g_lst; H_lst <- inte_data$H_lst
  theta_lst <- inte_data$theta_lst; n_lst <- inte_data$n_lst
  share_var <- inte_data$share_var; penalty.factor = inte_data$penalty
  
  XY_data <- pre_inte_XY(H_lst, g_lst, theta_lst, n_lst, 
                             share_var, penalty.factor)
  X_all <- XY_data$X_all
  Y_all <- XY_data$Y_all
  
  penalty.factor.all <- c(penalty.factor[which(share_var == 1)],
                          rep(penalty.factor[which(share_var == 0)], M))

  GIC.min <- Inf
  theta.opt <- NULL
  lambda.opt <- NULL
  
  for (lambda1 in lambda1_lst) {
    for (lambda2 in lambda2_lst) {
      penalty.factor.lambda <- penalty.factor.all
      penalty.factor.lambda[1:length(which(share_var == 1))] <- 
        penalty.factor.lambda[1:length(which(share_var == 1))] * lambda2
      fit.result <- integrative_fit_for_cv(X_all, Y_all, M, share_var, penalty.factor.lambda, 
                                           lambda1 = lambda1, lambda_init_beta = 0, gamma = gamma,
                                           max.iter = max.iter, tol = tol)
      GIC.lambda <- Cal_GIC(H_lst, g_lst, theta_lst, n_lst, share_var, fit.result, penalty.factor, 
                            lambda1 = lambda1, gamma = gamma, type = type)
      
      print(paste0("lam1=",round(lambda1,2)," lam2=",round(lambda2,2)))
      print(paste0("Value=",GIC.lambda))
      print("")
      if (GIC.lambda < GIC.min){
        GIC.min <- GIC.lambda
        theta.opt <- fit.result
        lambda.opt <- c(lambda1, lambda2, gamma)
      }
    }
  }
  print(paste0("Final: ", lambda.opt))
  return(list(theta = theta.opt, lambda.select = lambda.opt))
}
