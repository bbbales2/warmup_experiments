getHessian = function(fit, q) {
  Aqx = function(fit, q, r) {
    dx = 1e-5
    dr = dx * r
    (grad_log_prob(fit, q + dr / 2, adjust_transform = FALSE) -
        grad_log_prob(fit, q - dr / 2, adjust_transform = FALSE)) / dx
  }
  
  N = length(q)
  A = matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    r = rep(0, N)
    r[i] = 1
    A[, i] = Aqx(fit, q, r)
  }
  0.5 * (A + t(A))
}

diag_metric = function(samples) {
  diag(diag(cov(samples)))
}

dense_metric = function(samples) {
  c = cov(samples)
  
  if(nrow(samples) < ncol(samples)) {
    e = eigen(c, T)
    nkeep = nrow(samples) - 1
    mine = e$values[nkeep]
    c = e$vectors[, 1:nkeep] %*% diag(e$values[1:nkeep] - mine) %*% t(e$vectors[, 1:nkeep])
    c = c + mine * diag(ncol(samples))
  }
  
  return(c)
}

lw_linear_metric = function(samples) {
  linshrink_cov(samples)
}

lw_linear_corr_metric = function(samples) {
  sqrt_D = diag(sqrt(diag(cov(samples))))
  sqrt_Dinv = diag(1 / diag(sqrt_D))
  sqrt_D %*% linshrink_cov(samples %*% sqrt_Dinv) %*% sqrt_D
}

lw_nonlinear = function(X) {
  n = nrow(X)
  p = ncol(X)
  
  if(n < 12) {
    stop("n must be greater than 12")
  }
  
  sample = t(X) %*% X / n
  
  e = eigen(sample, T)
  u = e$vectors
  lambda = e$values
  
  isort = order(lambda)
  lambda = lambda[isort]
  u = u[, isort]
  
  lambda = lambda[max(1, p - n + 1) : p]
  L = matrix(rep(lambda, min(p, n)), nrow = length(lambda))
  h = n^(-1 / 3)
  H = h * t(L)
  
  x = (L - t(L)) / H
  
  ftilde = (3 / (4 * sqrt(5))) * rowMeans(pmax(1 - x^2 / 5, 0) / H)
  
  Hftemp = (-3 / (10 * pi)) * x +
    (3 / (4 * sqrt(5) * pi)) * (1 - x^2 / 5) * log(abs((sqrt(5) - x) / (sqrt(5) + x)))
  
  whichx = abs(x - sqrt(5)) <= 1e-16
  
  Hftemp[whichx] = (-3 / (10 * pi)) * x[whichx]
  Hftilde = rowMeans(Hftemp / H)
  
  if(p <= n) {
    dtilde = lambda / ((pi * (p / n) * lambda * ftilde)^2 +
                         (1 - (p / n) - pi * (p / n) * lambda * Hftilde)^2)
  } else {
    Hftilde0 = (1 / pi) * (3/(10 * h^2) + 3/(4 * sqrt(5) * h) * (1 - 1/(5 * h^2)) *
                             log((1 + sqrt(5) * h) / (1 - sqrt(5) * h))) * mean(1 / lambda)
    dtilde0 = 1 / (pi * (p - n) / n * Hftilde0)
    dtilde1 = lambda / (pi^2 * lambda^2 * (ftilde^2 + Hftilde^2))
    dtilde = c(rep(dtilde0, p - n), dtilde1)
  }
  
  u %*% diag(dtilde) %*% t(u)
}

lw_nonlinear_metric = function(samples) {
  #capture.output(out <- nlshrink_cov(samples))
  
  out = NULL
  
  tryCatch({
    out = lw_nonlinear(samples)
  }, error = function(e) {
    out = lw_linear_metric(samples)
    print(e)
  })
  
  return(out)
}

lw_nonlinear_corr_metric = function(samples) {
  sqrt_D = diag(sqrt(diag(cov(samples))))
  sqrt_Dinv = diag(1 / diag(sqrt_D))
  return(sqrt_D %*% lw_nonlinear_metric(samples %*% sqrt_Dinv) %*% sqrt_D)
}

lw_nonlinear2_metric = function(samples) {
  capture.output(out <- nlshrink_cov(samples))
  
  return(out)
}

lw_nonlinear2_corr_metric = function(samples) {
  sqrt_D = diag(sqrt(diag(cov(samples))))
  sqrt_Dinv = diag(1 / diag(sqrt_D))
  return(sqrt_D %*% lw_nonlinear2_metric(samples %*% sqrt_Dinv) %*% sqrt_D)
}

hessian_metric = function(samples) {
  H = getHessian(stan_fit, samples[1,])
  e = eigen(H, T)
  
  smallest_negative_evalue = max(e$values[e$values < 0])
  
  e$values[e$values >= 0] = smallest_negative_evalue
  
  return(e$vectors %*% diag(-1 / e$values) %*% t(e$vectors))
}

metrics = list(diag_metric = diag_metric,
               dense_metric = dense_metric,
               lw_linear_metric = lw_linear_metric,
               lw_linear_corr_metric = lw_linear_corr_metric,
               #lw_nonlinear_metric = lw_nonlinear_metric,
               #lw_nonlinear_corr_metric = lw_nonlinear_corr_metric,
               #lw_nonlinear2_metric = lw_nonlinear2_metric,
               #lw_nonlinear2_corr_metric = lw_nonlinear2_corr_metric,
               hessian_metric = hessian_metric)