com_psi <- function(t, knot){
  output <- cbind(1, t, matrix(0, length(t), length(knot) - 2))
  colnames(output) <- paste0('psi', 1:ncol(output))
  for (j in 2:(length(knot) - 1)){
    output[, paste0('psi', j + 1)] <-
      sapply(t, function(x){ x * ifelse(x > knot[j], 1, 0) })
  }
  return(t(output))
}

log_lkh <- function(eta, xi, beta, dat.obs, psi, d_phi, knot){
  ### compute full log likelihood of cox model
  cbhhat.obj <- cbhhat(dat.obs$y, eta, knot)
  phi <- t(dplyr::select(dat.obs, 'phi1':paste0('phi',d_phi)))
  nonx_col <- c('y', 'failed', 'calendar_time',
                paste0('phi', 1:d_phi))
  X <- t(dplyr::select(dat.obs, -nonx_col))
  eta <- matrix(eta, nrow = length(eta), ncol = 1)
  lglkh_1 <- dat.obs$failed *
    ( t(eta) %*% psi + t(xi) %*% phi + t(beta) %*% X )
  lglkh_2 <- -cbhhat.obj$int *
    exp( t(xi) %*% phi + t(beta) %*% X )
  lglkh <- -mean(lglkh_1 + lglkh_2)

  return(lglkh)
}




integral_target <- function(t_vec, knot, eta){
  psi <- com_psi(t_vec, knot)
  return(exp(t(eta %*% psi)))
}

cbhhat <- function(t_vec, eta, knot, n.nodes=15){
  gaussquad <- statmod::gauss.quad(n.nodes, kind = 'legendre')
  cbhazard <- sapply(1:n.nodes, function(i){
    result <- t_vec/2 * gaussquad$weights[i] *
      integral_target(t_vec * (gaussquad$nodes[i] + 1) / 2, knot, eta)
    return(result)
  })

  int <- rowSums(cbhazard)
  return(list(int = int))
}



grr <- function(eta, xi, beta, dat.obs, psi, d_phi, knot, n.nodes = 15){
  ### auxiliary function for optim, which compute the gradient of our -log-likelihood
  t_vec = dat.obs$y
  gaussquad <- statmod::gauss.quad(n.nodes, kind = 'legendre')
  phi <- t(dplyr::select(dat.obs, 'phi1':paste0('phi', d_phi)))
  nonx_col <- c('y','failed','calendar_time',paste0('phi', 1:d_phi))
  X <- t(dplyr::select(dat.obs, -nonx_col))
  Delta <- dat.obs$failed
  Z <- as.numeric(t(xi) %*% phi + t(beta) %*% X)

  grr_part <- sapply(1:n.nodes, function(i){
    result <- t_vec/2 * gaussquad$weights[i] *
      as.vector(integral_target(t_vec * (gaussquad$nodes[i] + 1) / 2, knot, eta)) *
      t(com_psi(t_vec * (gaussquad$nodes[i] + 1) / 2, knot))
    return(list(result))
  })

  grr_part_i <- Reduce("+", grr_part)
  grad <- -1 * colSums((Delta * t(psi)) - exp(Z) * grr_part_i) / nrow(dat.obs)
  return(grad)
}

