global_IPD_loss <- function(fit.result, local_data_lst, gamma, lambda1, lambda2){
  loss.result = sum(sapply(1:dim(fit.result)[2], function(m){
    nrow(local_lst[[m]]$dat.obs) *
      log_lkh(eta = fit.result[str_detect(rownames(fit.result),"^psi"),m],
              xi =  as.matrix(fit.result[str_detect(rownames(fit.result),"^phi[0-9]+$"),m],
                              ncol = 1),
              beta =  as.matrix(fit.result[!str_detect(rownames(fit.result),"^phi[0-9]+$") &
                                             !str_detect(rownames(fit.result),"^psi"),m],
                                ncol = 1),
              dat.obs = local_lst[[m]]$dat.obs,
              psi = local_lst[[m]]$psi,
              d_phi = sum(str_detect(rownames(fit.result),"^phi[0-9]$")),
              knot = local_lst[[m]]$psi.knot
      ) 
  }))
  
  pen1 = fit.result[str_detect(rownames(fit.result),"^phi[0-9]+\\."), ]
  pen1.sub.sum = c()
  for(phi.j in 1:sum(str_detect(rownames(fit.result),"^phi[0-9]$"))){
    pen1.sub = pen1[str_detect(rownames(pen1),as.character(phi.j)),]
    pen1.sub.sum = cbind(pen1.sub.sum,
                         as.matrix(apply(pen1.sub,1,function(x){
                           sum(abs(x))
                         })))
  }
  pen1.sum = sum(apply(pen1.sub.sum, 1, function(x){
    sum(abs(x))^gamma
  }))   ## pen1 for phi
  
  pen2 = fit.result[!str_detect(rownames(fit.result),"^phi") &
                      !str_detect(rownames(fit.result),"^psi"),1]
  loss.result = loss.result + lambda2 * lambda1 * sum(abs(pen2)) + lambda1 * pen1.sum
  return(loss.result)
}



eval_fun <- function(dat_obj, local_est, theta_fit){
  # dat_obj is a testing data object containing the testing samples
  # and their true risk given the features and calendar times.
  
  # local_est is a processed data of dat_obj with the spline basis and interaction features
  # (so that we can use the fitted (whole)coefficients vector theta_fit)
  
  # Calculate the prediciton:
  
  X_test <- cbind(t(local_est$psi), local_est$phi, t(local_est$X))
  pred_val <- as.vector(X_test %*% theta_fit)
  
  # Calculate the error of the predicted risk (lp, linear predictor since required by AUC.sh)
  
  pred_loss <- mean((dat_obj$lp - mean(dat_obj$lp) - pred_val + mean(pred_val))^2)
  
  # AUC based on the test data in dat_obj
  
  train_indx <- 1:as.integer(length(dat_obj$y) / 2)
  test_indx <- setdiff(1:length(dat_obj$y), train_indx)
  auc_eval <- AUC.sh(Surv.rsp = Surv(dat_obj$y[train_indx], dat_obj$failed[train_indx]), 
                     Surv.rsp.new = Surv(dat_obj$y[test_indx], dat_obj$failed[test_indx]),
                     lp = pred_val[train_indx], lpnew = pred_val[test_indx], 
                     times = seq(min(dat_obj$y), median(dat_obj$y), length.out = 10))
  
  return(list(pred_loss = pred_loss,
              #nllh = nlog_llh,
              iauc = auc_eval$iauc,
              auc_start = auc_eval$auc[1],
              auc_median = auc_eval$auc[5]
              #beta_error = beta_error
  )
  )
}
