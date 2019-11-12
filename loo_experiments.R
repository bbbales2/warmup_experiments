library(mvtnorm)
library(cmdstanr)
library(tidyverse)
library(rstan)
library(nlshrink)
library(parallel)

set_cmdstan_path("/home/bbales2/cmdstan-warmup/")

models = tibble(name = c("kilpisjarvi", "accel_splines", "accel_gp", "diamonds", "prophet", "radon"),
                adapt_delta = c(0.8, 0.99, 0.99, 0.8, 0.8, 0.8),
                init = c(2.0, 0.25, 0.25, 2.0, 2.0, 2.0))

out = list()
for(model_i in 1:nrow(models)) {
  model_name = models[[model_i, "name"]]
  
  model_path = paste0(cmdstan_path(), "/examples/", model_name, "/", model_name, ".stan")
  
  model = cmdstan_model(model_path)
  model$compile()
  
  data_env = new.env(parent = baseenv())
  source(paste0(cmdstan_path(), "/examples/", model_name, "/", model_name, ".dat"), local = data_env)
  data = as.list.environment(data_env)
  
  stan_fit = stan(model_path,
                  data = data,
                  iter = 1)
  
  fit = model$sample(data = data,
                     num_chains = 1,
                     save_warmup = 1,
                     experimental = 1,
                     metric = "dense_e",
                     which_adaptation = 0,
                     save_diagnostics = TRUE,
                     init = models[[model_i, "init"]],
                     adapt_delta = models[[model_i, "adapt_delta"]])
  
  ref_fit = model$sample(data = data,
                         num_chains = 4,
                         num_samples = 2000,
                         save_warmup = 0,
                         experimental = 1,
                         metric = "dense_e",
                         which_adaptation = 0,
                         init = models[[model_i, "init"]],
                         adapt_delta = models[[model_i, "adapt_delta"]],
                         save_diagnostics = TRUE)
  
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
  getUnconstrainedSamples = function(files) {
    lapply(files, function(file) {
      read.csv(file, comment.char = "#") %>%
        as_tibble() %>%
        select(-lp__, -accept_stat__, -stepsize__, -treedepth__, -n_leapfrog__, -divergent__, -energy__,
               -starts_with("p_"), -starts_with("g_"))
    }) %>%
      bind_rows() %>%
      as.matrix()
  }
  
  unconstrained_samples = getUnconstrainedSamples(fit$diagnostic_files())
  ref_unconstrained_samples = getUnconstrainedSamples(ref_fit$diagnostic_files())
  ref_C = cov(ref_unconstrained_samples)
  
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
  
  metric_df = lapply(c(5, 10, 15, 20, 25, 30, 35, 40, 80), function(window_length) {#, 160, 320
    print(paste0("window_length: ", window_length))
    
    lapply(names(metrics), function(metric_name) {
      lpds = sapply(1:16, function(i) {
        istart = sample(1:(nrow(unconstrained_samples) - window_length), 1)
        iend = istart + window_length
        
        y = unconstrained_samples[istart:iend, ]
        
        tryCatch({
          lpds = c()
          for(loi in 1:window_length) {
            yp = y[loi, ]
            yloo = y[-loi, ]
            metric = metrics[[metric_name]](yloo)
            lpds = c(lpds, dmvnorm(yp, colMeans(yloo), metric, log = TRUE))
          }
        }, error = function(e) {
          print(e)
          lpds = NA
        })
        
        return(mean(lpds))
      })
      
      print(paste0(model_name, ", metric_name: ", metric_name, ", window_length: ", window_length))
      
      return(tibble(model_name = model_name,
                    lpd = lpds,
                    metric = metric_name,
                    window_length = window_length))
    }) %>% bind_rows
  }) %>% bind_rows
  
  out[[model_i]] = metric_df
}

out %>%
  bind_rows() %>%
  filter(lpd > -200) %>%
  ggplot() +
  geom_jitter(aes(window_length, lpd, color = metric), size = 1.0, alpha = 1.0, width = 0.1) +
  scale_x_log10() +
  #scale_y_log10() +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1.0))) +
  facet_grid(~ model_name)

istart = sample(1:(nrow(unconstrained_samples) - window_length), 1)
iend = istart + window_length

y = unconstrained_samples[istart:iend, ]

for(metric_name in names(metrics)) {
  lpds = c()
  for(loi in 1:window_length) {
    yp = unconstrained_samples[iend + 10 + loi, ]
    #yp = y[loi, ]
    yloo = y[-loi, ]
    metric = metrics[[metric_name]](yloo)
    if(!is.null(metric)) {
      #mu = colMeans(yloo)
      #lpds = c(lpds, -(yp - mu) %*% solve(metric, yp - mu) - log(det(metric)) - length(mu) * log(2 * pi))
      lpds = c(lpds, dmvnorm(yp, colMeans(yloo), metric, log = TRUE))
    }
  }
  
  metric = metrics[[metric_name]](yloo)
  L = t(chol(metric))
  e = eigen(solve(L, t(solve(L, t(ref_C)))), T)
  cn = max(abs(e$values)) / min(abs(e$values))
  
  stan_rdump("inv_metric", paste0(metric_name, ".dat"), env = list2env(list(inv_metric = metric)))
  
  print(paste(metric_name, mean(lpds), cn, sd(lpds) / sqrt(length(lpds))))
}

fit2 = model$sample(data = data,
                   num_chains = 1,
                   save_warmup = 1,
                   experimental = 1,
                   metric = "dense_e",
                   inv_metric = dense_metric(unconstrained_samples),
                   which_adaptation = 0,
                   save_diagnostics = TRUE,
                   init = models[[model_i, "init"]],
                   adapt_delta = models[[model_i, "adapt_delta"]])

fit2 = model$sample(data = data,
                   num_chains = 1,
                   save_warmup = 0,
                   metric = "diag_e",
                   inv_metric = diag(diag_metric(y)),
                   init = 0,
                   adapt_engaged = 1,
                   num_warmup = 250,
                   init_buffer = 0,
                   term_buffer = 250,
                   window = 0)

fit2$summary()

ldf = lapply(names(metrics), function(metric_name) {
  metric = metrics[[metric_name]](y)
  rands = rmvnorm(1000, colMeans(y), metric)
  tibble(lp = apply(rands, 1, function(x) log_prob(stan_fit, x)),
         lq = dmvnorm(rands, colMeans(y), metric, log = TRUE),
         name = metric_name)
}) %>% bind_rows

ldf %>%
  mutate(lr = lp - lq) %>%
  group_by(name) %>%
  summarize(m = median(lr),
            l = quantile(lr, 0.2),
            u = quantile(lr, 0.8))

pdf = lapply(names(metrics), function(metric_name) {
  fit = model$sample(data = data,
                     num_chains = 4,
                     metric = "dense_e",
                     inv_metric = metrics[[metric_name]](y),
                     init = 0,
                     adapt_engaged = 1,
                     num_warmup = 250,
                     init_buffer = 0,
                     term_buffer = 250,
                     window = 0)
  
  process_fit = function(files, draws) {
    stanfit = rstan::read_stan_csv(files)
    sp = rstan::get_sampler_params(stanfit)
    fitstats = monitor(draws, warmup = 0, print = FALSE)
    
    efficiency_df = tibble(rhat = max(fitstats$Rhat),
                           min_bulk_ess = min(fitstats$Bulk_ESS),
                           min_tail_ess = min(fitstats$Tail_ESS),
                           lp_bulk_ess = fitstats[["lp__", "Bulk_ESS"]],
                           lp_tail_ess = fitstats[["lp__", "Tail_ESS"]],
                           nleapfrogs = sapply(sp, function(e) {
                             e[, "n_leapfrog__"]
                           }) %>% sum,
                           ndivergences = sapply(sp, function(e) {
                             e[, "divergent__"]
                           }) %>% sum)
    
    return(efficiency_df)
  }
  
  process_fit(fit$output_files(), fit$draws()) %>%
    mutate(method = metric_name)
}) %>%
  bind_rows()

pdf %>%
  mutate(name = names(metrics)) %>%
  mutate(min_bulk_ess / nleapfrogs)

read_csv(fit2$output_files(), comment = "#") %>%
  mutate(r = row_number()) %>%
  filter(r < 250) %>%
  ggplot(aes(r, lp__)) +
  geom_point()

dense_metric(y)


model = cmdstan_model("beta.stan")
fit = model$sample(data = data,
                   num_chains = 1,
                   save_warmup = 1,
                   save_diagnostics = TRUE)


pm = read_csv(fit$diagnostic_files(), comment = "#") %>%
  select(starts_with("p_")) %>%
  as.matrix()

diagf = read_csv(fit$diagnostic_files(), comment = "#")

gm = read_csv(fit$diagnostic_files(), comment = "#") %>%
  select(starts_with("g_")) %>%
  as.matrix()

gm %*% t(L)

tibble(v = sapply(1:nrow(pm), function(i) {
  v = acos(abs((pm[i, ] %*% gm[i, ]) / (sqrt(sum(pm[i, ]^2)) * sqrt(sum(gm[i, ]^2))))) * 180 / pi
  return(v)
})) %>%
  mutate(r = row_number()) %>%
  ggplot(aes(r, v)) +
  geom_point()

tibble(sum_gradients = rowSums(gm)) %>%
  mutate(r = row_number(),
         lp = diagf$lp__) %>%
  filter(r < 250) %>%
  gather(which, value, -r) %>%
  ggplot(aes(r, value)) +
  geom_point(aes(color = which))


gm2 = read_csv(fit2$diagnostic_files(), comment = "#") %>%
  select(starts_with("g_")) %>%
  as.matrix()

gm3 = t(apply(unconstrained_samples, 1, function(row) grad_log_prob(stan_fit, row)))

lb = sqrt(qchisq(0.025, ncol(gm)))
ub = sqrt(qchisq(0.975, ncol(gm)))

sqrtp = sqrt(rowSums(pm[1:100,]^2))

qbinom(0.05, 8, 0.95) <
  ((lb < sqrtp & sqrtp < ub) %>%
         rollapply(8, sum))

tibble(magnitude_momentum = sqrt(rowSums(pm^2)),
       sum_gradients = rowSums(gm)) %>%
  mutate(draw = row_number(),
         lp = read_csv(fit$diagnostic_files(), comment = "#")$lp__) %>%
  filter(draw < 100 & lp > -500) %>%
  gather(which, value, -draw) %>%
  group_by(which) %>%
  mutate(unitless_scale = value / sd(value)) %>%
  ggplot(aes(draw, unitless_scale)) +
  geom_point(aes(color = which))

tibble(v1 = rowSums(gm %*% diag(sqrt(diag(diag_metric(unconstrained_samples))))),
       v2 = rowSums(gm %*% t(L))) %>%
  mutate(r = row_number(),
         lp = diagf$lp__) %>%
  filter(r < 100 & lp > -500) %>%
  gather(which, value, -r) %>%
  ggplot(aes(r, value)) +
  geom_point(aes(color = which))

diagf %>%
  mutate(draw = row_number()) %>%
  filter(draw < 100 & energy__ < 200) %>%
  ggplot(aes(draw, n_leapfrog__)) +
  geom_point(aes(color = as.factor(divergent__)))

diagf %>%
  select(starts_with("g_")) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = c("variable"),
                  names_pattern = "g_(.*)",
                  values_to = "g") %>%
  filter(draw < 100) %>%
  filter(abs(g) < 10000) %>%
  ggplot(aes(draw, g)) +
  geom_point(aes(color = variable), alpha = 0.5, size = 0.5) +
  theme(legend.position = "none")

tibble(v = sapply(1:nrow(pm), function(i) {
  v = sum(gm[i,1])
  return(v)
}))

mean(gm[1000:2000,])

hist(gm[1000:2000])

gm[1500:2000,82] %>% sd

%>% mean
            
            %>%
  mutate(r = row_number()) %>%
  pivot_longer(-r, names_to = c("variable", "idx"),
               names_pattern = "p_(.*).([0-9]+)",
               values_to = "v")

