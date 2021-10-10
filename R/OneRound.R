Eval_local_inte_IPD <- function(seed){

  ################################################################################
  ##########################  Local data generation ##############################
  ################################################################################

  K = 5
  n_vec = c(3000,2000,1000,1500,2000)
  T_vec = c(300,210,115,150,250)
  beta_nonzero_num = 6
  beta_zero_num = 44

  local_raw_lst = generate_data_Multilocals(seed = 123, n_vec, T_vec,
                                            beta_nonzero_num, beta_zero_num)


  ################################################################################
  ##########################  Local par estimate #################################
  ################################################################################

  t.local.start <- proc.time()
  local_lst = c()
  for(k in 1:K){
    print(paste0("Site ", k))
    #### local fit:
    local <- surv_local(local_raw_lst[[k]],
                        num_bh_knot=2, d_phi=3, n.nodes=15)
    sum_gH <- summary_gH(local$dat.obs, local$local_est)
    sum_stat <- c(`dat.obs` = list(local$dat.obs), local$local_est, sum_gH)
    proc.time() - t
    local_lst <- c(local_lst, list(sum_stat))
  }
  t.local <- proc.time() - t.local.start


  ################################################################################
  ##########################  Integrative par estimate ###########################
  ################################################################################

  t.inte.start <- proc.time()
  inte_data <- pre_inte_fromLocal(local_lst)

  # Hyper par
  lambda1_lst <- sqrt(mean(inte_data$n_lst) * log(length(inte_data$g_lst[[1]]))) *
    exp(seq(from = log(0.05), to = log(0.5), length.out = 5))
  lambda2_lst <- exp(seq(from = log(1), to = log(30), length.out = 5))   # LASSO
  gamma <- 0.5


  cv.fit.aggregate <- integrative_fit(inte_data,
                                      lambda1_lst = lambda1_lst, lambda2_lst = lambda2_lst,
                                      gamma = gamma, type = 'BIC')

  t.inte <- proc.time() - t.inte.start

  ################################################################################
  #############################  IPD par estimate #################################
  ################################################################################

  t.IPD.start <- proc.time()
  cv.IPD.aggregate <- IPD_fit(inte_data, local_lst,
                              lambda1_lst = lambda1_lst, lambda2_lst = lambda2_lst,
                              gamma = gamma, type = 'BIC')

  t.IPD <- proc.time() - t.IPD.start

  ################################################################################
  ##################################  Evaluation #################################
  ################################################################################

  eval_results <- vector('list', K)
  for (k in 1:K){
    dat_obj <- local_raw_lst[[k]]
    local_est <- local_lst[[k]]
    theta_loc <- c(as.vector(local_est$eta), as.vector(local_est$xi), as.vector(local_est$beta))
    theta_int <- cv.fit.aggregate$theta[,k]
    theta_IPD <- cv.IPD.aggregate$theta[,k]

    eval_num <- list(local = eval_fun(dat_obj, local_est, theta_loc),
                     our = eval_fun(dat_obj, local_est, theta_int),
                     ipd = eval_fun(dat_obj, local_est, theta_IPD))

    eval_results[[k]] <- eval_num
  }
  return(eval_results)
}

