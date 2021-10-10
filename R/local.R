#==============================================================================================
#=================================== Local parameter estimation ===============================
#==============================================================================================

surv_local <- function(local_raw, num_bh_knot=2, d_phi=3, n.nodes=15){
  ### the main function for local parameter estimation and summary statistics derivation
  ### input: time: in numeric, survival time;
  ###        x: in numeric, clinical feature (if value is 0 and included in missing_col,
  ###            will construct missing indicator for the corresponding column) ;
  ###        event: in binary, censorship (1 for not censored);
  ###        calendar_time: in numeric, calendar_time, used to construct phi
  ###        num_bh_knot: in int, number of nodes of piecewise-exponental function we use to
  ###                     approximate baseline hazard, num_bh_knot=dimension of psi, we can
  ###                     simply use the default value 2 nodes, corresponding to dim(psi)=4
  ###        d_phi: in tin, number of degrees of freedom of cubic spline for calendar time,
  ###                we can use default value, 3
  ### output: a list containing estimations for all parameters and summary statistics,
  ###                 including gradient and Hessian

  x <- local_raw[, stringr::str_detect(colnames(local_raw), "X")]
  time <- local_raw$y
  event <- local_raw$failed
  calendar_time <- local_raw$calendar_time

  # add cubic baseline for calendar time
  phi <- data.frame(splines::ns(calendar_time, df = d_phi))
  colnames(phi) <- paste0('phi', 1:ncol(phi))

  # make it shifted to be orthogonal to 1
  for (j in 1:d_phi) {
    phi[, j] <- (phi[, j] - mean(phi[, j])) / stats::sd(phi[, j])
  }

  # interaction with covariates
  X_ob <- x
  for (j in 1:d_phi) {
    inter_basis <- phi[,j] * x
    colnames(inter_basis) <- paste0(colnames(phi)[j], ".", colnames(x))
    X_ob <- X_ob %>% dplyr::mutate(inter_basis)
  }
  # interaction with missing indicator
  dat.obs <- dplyr::mutate(cbind(X_ob, phi),
                    y = time,
                    failed = event,
                    calendar_time = calendar_time)
  local_est <- get_local_est(dat.obs, num_bh_knot, d_phi, n.nodes)
  return(list(`dat.obs` = dat.obs,
              `local_est` = local_est))
}



get_local_est <- function(dat.obs, num_bh_knot, d_phi, n.nodes = 15){
  ### get estimation of parameters at local sites
  ### input: dat.obs: the covariate matrix;
  ###        num_bh_knot: number of nodes used in piecewise-exponential function to approximate the baseline hazard function
  ###        d_phi: dimensions for phi, the degrees of freedom of the cubic splines of calendar time
  ### output: all estimations needed to compute summary statistics

  nonx_col <- c('y', 'failed', 'calendar_time')
  xdata <- dplyr::select(dat.obs, -dplyr::all_of(nonx_col))

  # fit stratified adaptive lasso
  penalty0 <- data.frame(matrix(1, 1, ncol(xdata)))
  colnames(penalty0) <- colnames(xdata)
  penalty0[paste0('phi',1:d_phi)] = 0  # no penalty for calendar effect

  # first fit a ridge regression
  ridge_model <- glmnet::cv.glmnet(x = as.matrix(xdata),
                          y = survival::Surv(dat.obs$y, dat.obs$failed),
                          family = 'cox',
                          alpha = 0,
                          penalty.factor = penalty0)
  # adaptivelasso, use absolute value coefficients of ridge as reciprocal of lasso penalty
  penalty1 = as.vector(penalty0) / as.vector(abs(stats::coef(ridge_model,s='lambda.min')))
  alasso_model = glmnet::cv.glmnet(x = as.matrix(xdata),
                           y = survival::Surv(dat.obs$y, dat.obs$failed),
                           family = 'cox',
                           alpha = 1,
                           penalty.factor = as.vector(penalty1))
  mymodel <- stats::as.formula(paste0('Surv(y, failed) ~',
                               paste0(colnames(xdata), collapse='+')))
  cox_model <- survival::coxph(mymodel, data = dat.obs)
  cox_model$coefficients <- stats::coef(alasso_model, s = 'lambda.min')

  # extract beta
  coeff_theta <- stats::coef(alasso_model, s = 'lambda.min')
  betahat <- as.matrix(coeff_theta[setdiff(colnames(xdata), c(paste0('phi', 1:d_phi))), ])
  phi <- as.matrix(dplyr::select(dat.obs, paste0('phi', 1:d_phi)))
  xi <- as.matrix(coeff_theta[paste0('phi', 1:d_phi), ])

  # estimate eta for baseline hazard
  psi.knot <- stats::quantile(dat.obs$y, seq(0, 1, length.out = num_bh_knot + 2))
  psi.knot[1] <- 0
  psi.knot[length(psi.knot)] <- Inf
  psi <- com_psi(dat.obs$y, psi.knot)

  eta_initial <- rep(.1, nrow(psi))
  eta <- stats::optim(eta_initial,
               fn = log_lkh, xi = xi, beta = betahat, dat.obs = dat.obs,
               psi = psi, d_phi = d_phi, knot = psi.knot,
               gr = grr, method = 'L-BFGS-B', lower = c(rep(-Inf, length(eta_initial)))
               )$par


  cbhhat.obj <- cbhhat(dat.obs$y, eta, psi.knot, n.nodes)
  myX <- as.matrix(dplyr::select(xdata, -paste0('phi', 1:d_phi)))
  return(list(X = t(myX),
              beta = betahat,
              eta = eta,
              xi = xi,
              phi = phi,
              psi = psi,
              psi.knot =  psi.knot,
              n.nodes = n.nodes,
              alasso = alasso_model,
              cbhhat.obj = cbhhat.obj))
}



summary_gH <- function(dat.obs, local_est){
  eta <- local_est$eta
  xi <- local_est$xi
  beta <- local_est$beta
  psi <- local_est$psi
  knot <- local_est$psi.knot
  d_phi <- ncol(local_est$phi)
  nn <- length(dat.obs$y)
  t_vec <- dat.obs$y
  n.nodes <- local_est$n.nodes
  gaussquad <- statmod::gauss.quad(n.nodes, kind = 'legendre')
  phi <- t(dplyr::select(dat.obs, 'phi1':paste0('phi', d_phi)))
  nonx_col <- c('y','failed','calendar_time',paste0('phi', 1:d_phi))
  X <- t(dplyr::select(dat.obs, -nonx_col))
  Delta <- dat.obs$failed
  beta_full_Z <- as.numeric(t(xi) %*% phi + t(beta) %*% X)
  Z <- t(rbind(phi, X))
  g_eta_part = g_beta_part =
    H_eta_eta_part = H_eta_beta_part = H_beta_beta_part = 0
  for(i in 1:n.nodes){
    phi_t_new_q <- t(com_psi(t_vec * (gaussquad$nodes[i] + 1) / 2, knot))
    f_i_q <- t_vec/2 * gaussquad$weights[i] *
      exp(matrix(eta, 1, length(eta)) %*% t(phi_t_new_q) + beta_full_Z)
    g_eta_part_i <- sweep(phi_t_new_q, as.vector(f_i_q), MARGIN = 1, '*')
    g_beta_part_i <- sweep(Z, as.vector(f_i_q), MARGIN = 1, '*')
    H_eta_eta_part <- H_eta_eta_part + t(g_eta_part_i) %*% phi_t_new_q
    H_eta_beta_part <- H_eta_beta_part + t(g_eta_part_i) %*% Z
    H_beta_beta_part <- H_beta_beta_part + t(g_beta_part_i) %*% Z
    g_eta_part = g_eta_part + colSums(g_eta_part_i)
    g_beta_part = g_beta_part + colSums(g_beta_part_i)
  }
  g_eta <- g_eta_part / nn -
    colMeans(sweep(t(com_psi(t_vec, knot = knot)), Delta, MARGIN = 1, '*'))
  g_beta <- g_beta_part / nn -
    colMeans(sweep(Z, Delta, MARGIN = 1, '*'))
  H_eta_eta <- H_eta_eta_part / nn
  H_eta_beta <- H_eta_beta_part / nn
  H_beta_beta <- H_beta_beta_part / nn
  grad <- c(g_eta, g_beta)
  H <- rbind(cbind(H_eta_eta, H_eta_beta),
             cbind(t(H_eta_beta), H_beta_beta))
  return(list(`grad` = grad,
              `Hessian` = H))
}


