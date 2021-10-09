IPD_fit <-  function(inte_data,
                        local_lst,
                        lambda1_lst = 0.1, lambda2_lst = lambda1_lst,
                        lambda_init_beta = 0, gamma = 0.5,
                        max.iter = 10, tol = 1e-5,
                        type = 'BIC'){
  ## Take local data
  M <- length(inte_data$n_lst)
  g_lst <- inte_data$g_lst; H_lst <- inte_data$H_lst
  theta_lst <- inte_data$theta_lst; n_lst <- inte_data$n_lst
  share_var <- inte_data$share_var; penalty.factor = inte_data$penalty

  GIC.min <- Inf
  theta.opt <- NULL
  lambda.opt <- NULL
  penalty.factor.all <- c(penalty.factor[which(share_var == 1)],
                          rep(penalty.factor[which(share_var == 0)], M))
  XY_data <- pre_inte_XY(H_lst, g_lst, theta_lst, n_lst, 
                         share_var, penalty.factor)
  X_all0 <- XY_data$X_all
  Y_all0 <- XY_data$Y_all
  
  for (lambda1 in lambda1_lst) {
    for (lambda2 in lambda2_lst) {
      print(paste0("lambda1: ", lambda1))
      print(paste0("lambda2: ", lambda2))
      penalty.factor.lambda <- penalty.factor.all
      penalty.factor.lambda[1:length(which(share_var == 1))] <- penalty.factor.lambda[1:length(which(share_var == 1))] * lambda2
      
      fit.result.new <- fit.result <- integrative_fit_for_cv(X_all0, Y_all0, M, share_var, penalty.factor.lambda, 
                                           lambda1 = lambda1, lambda_init_beta = 0, gamma = gamma,
                                           max.iter = max.iter, tol = tol)
      
      GIC.lambda.old <- Cal_GIC(H_lst, g_lst, theta_lst, 
                                n_lst, share_var, fit.result, penalty.factor, 
                                lambda1 = lambda1, gamma = gamma, type = type)
      print(paste0("integrative GIC: ", GIC.lambda.old))
      
      loss.result = global_IPD_loss(fit.result, local_lst, gamma, lambda1, lambda2)
      ii = 0
      loss.change = -1
      local_lst_new = local_lst
      H_lst.new <- H_lst
      g_lst.new <- g_lst
      theta_lst.new <- theta_lst
      
      while(loss.change < 0){
        ii = ii + 1
        print(ii)
        loss.prev = loss.result
        H_lst_lambda <- H_lst.new
        g_lst_lambda <- g_lst.new
        theta_lst_lambda <- theta_lst.new
        fit.result_lambda <- fit.result.new
        for(m in 1:M){
          local_lst_new[[m]]$eta = 
            fit.result.new[str_detect(rownames(fit.result.new),"^psi"),m]
          local_lst_new[[m]]$xi = 
            as.matrix(fit.result.new[str_detect(rownames(fit.result.new),"^phi[0-9]+$"),m],
                      ncol = 1)
          local_lst_new[[m]]$beta = 
            as.matrix(fit.result.new[!str_detect(rownames(fit.result.new),"^phi[0-9]+$") &
                                   !str_detect(rownames(fit.result.new),"^psi"),m],
                      ncol = 1)
          a = summary_gH(dat.obs = local_lst_new[[m]]$dat.obs, 
                         local_est = local_lst_new[[m]]) 
          local_lst_new[[m]]$grad = g_lst.new[[m]] = a$grad
          local_lst_new[[m]]$Hessian = H_lst.new[[m]] = a$Hessian
          theta_lst.new[[m]] <- c(local_lst_new[[m]]$eta,
                                  local_lst_new[[m]]$xi, 
                                  local_lst_new[[m]]$beta)
        }
        
        XY_data <- pre_inte_XY(H_lst.new, g_lst.new, theta_lst.new, n_lst, 
                               share_var, penalty.factor)
        X_all <- XY_data$X_all
        Y_all <- XY_data$Y_all
        
        fit.result.new <- integrative_fit_for_cv(X_all, Y_all, M, share_var, penalty.factor.lambda, 
                                             lambda1 = lambda1, lambda_init_beta = 0, gamma = gamma,
                                             max.iter = max.iter, tol = tol)

        loss.result = global_IPD_loss(fit.result.new, local_lst_new, gamma, lambda1, lambda2)
        loss.change = loss.result - loss.prev
        print(paste0("Loss change: ",loss.change))
        GIC.lambda.new <- Cal_GIC(H_lst.new, g_lst.new, theta_lst.new, 
                              n_lst, share_var, fit.result.new, penalty.factor, 
                              lambda1 = lambda1, gamma = gamma, type = type)
        print(paste0("Updated GIC: ", GIC.lambda.new))
        # fit.change = (norm(fit.result[!str_detect(rownames(fit.result),"^phi[0-9]+$") &
        #                                 !str_detect(rownames(fit.result),"^psi"),] - 
        #                      fit.result[!str_detect(rownames(fit.result),"^phi[0-9]+$") &
        #                      !str_detect(rownames(fit.result),"^psi"),], type = '2'))^2
        # 
        # 
        # fit.change2 = (norm(as.vector(fit.result) - as.vector(fit.result), type = '2'))^2
      }
      
      GIC.lambda <- Cal_GIC(H_lst_lambda, g_lst_lambda, theta_lst_lambda, 
                            n_lst, share_var, fit.result_lambda, penalty.factor, 
                            lambda1 = lambda1, gamma = gamma, type = type)
      print(paste0("GIC lambda: ",GIC.lambda))
      
      

      
      
      print("-------------------------------------")
      print("-------------------------------------")
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

