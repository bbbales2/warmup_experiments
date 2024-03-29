library(mvtnorm)
library(cmdstanr)
library(tidyverse)
library(rstan)
library(nlshrink)
library(parallel)

set_cmdstan_path("/home/bbales2/warmup/cmdstan")

models = tibble(name = c("kilpisjarvi", "accel_splines", "accel_gp", "diamonds", "prophet", "radon"),
                adapt_delta = c(0.8, 0.99, 0.99, 0.8, 0.8, 0.8),
                init = c(2.0, 0.25, 0.25, 2.0, 2.0, 2.0))

out = list()
#for(model_i in 1:nrow(models)) {
model_i = 4

model_name = models[[model_i, "name"]]
  
model_path = paste0(cmdstan_path(), "/examples/", model_name, "/", model_name, ".stan")
  
model = cmdstan_model(model_path, quiet = FALSE)

data_env = new.env(parent = baseenv())
source(paste0(cmdstan_path(), "/examples/", model_name, "/", model_name, ".dat"), local = data_env)
data = as.list.environment(data_env)
  
stan_fit = stan(model_path,
                data = data,
                iter = 1)

diag_inv_metric = function(samples) {
  diag(diag(cov(samples)))
}
dense_inv_metric = function(samples) {
  c = cov(samples)
  
  if(nrow(samples) < ncol(samples)) {
    e = eigen(c, T)
    nkeep = nrow(samples) - 1
    evalues = pmax(e$values, 1e-16)
    mine = evalues[nkeep]
    c = e$vectors[, 1:nkeep] %*% diag(evalues[1:nkeep] - mine) %*% t(e$vectors[, 1:nkeep])
    c = c + mine * diag(ncol(samples))
  }
  
  return(c)
}
lw_linear_inv_metric = function(samples) {
  linshrink_cov(samples)
}
lw_linear_corr_inv_metric = function(samples) {
  sqrt_D = diag(sqrt(diag(cov(samples))))
  sqrt_Dinv = diag(1 / diag(sqrt_D))
  sqrt_D %*% linshrink_cov(samples %*% sqrt_Dinv) %*% sqrt_D
}

inv_metrics = list(diag = diag_inv_metric,
                   dense = dense_inv_metric,
                   lw_linear = lw_linear_inv_metric,
                   lw_linear_corr = lw_linear_corr_inv_metric)

# Get in typical set
fit = model$sample(data = data,
                   num_chains = 1,
                   save_warmup = 1,
                   num_warmup = 100,
                   num_sample = 0,
                   experimental = 1,
                   metric = "dense_e",
                   which_adaptation = 0,
                   save_diagnostics = TRUE,
                   init = models[[model_i, "init"]],
                   adapt_delta = models[[model_i, "adapt_delta"]],
                   init_buffer = 100,
                   window = 0,
                   term_buffer = 0)

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

getLeapfrogs = function(files) {
  sapply(files, function(file) {
    read.csv(file, comment.char = "#") %>%
      as_tibble() %>%
      pull(n_leapfrog__) %>%
      sum()
  }) %>%
    sum
}

getExtras = function(files) {
  lapply(files, function(file) {
    read.csv(file, comment.char = "#") %>%
      as_tibble() %>%
      select(lp__, accept_stat__, stepsize__, treedepth__, n_leapfrog__, divergent__, energy__)
  }) %>%
    bind_rows()
}

usamples = getUnconstrainedSamples(fit$diagnostic_files())

get_init = function(usamples) {
  ldraw = usamples %>%
    tail(1)

  init = constrain_pars(stan_fit, ldraw %>% as.matrix)
  init_file = tempfile("init", fileext = ".dat")
  stan_rdump(names(init), init_file, env = list2env(init))
  
  return(init_file)
}

get_stepsize = function(fit) {
  read_csv(fit$diagnostic_files(), comment = "#") %>%
    tail(1) %>%
    pull(stepsize__)
}

inv_metric = diag(ncol(usamples))

Nw = 50
#get_slow_adapting_samples = function(Nw, q, inv_metric) {

init_file = get_init(usamples)
stepsize = get_stepsize(fit)
fit = model$sample(data = data,
                   num_chains = 1,
                   save_warmup = 1,
                   num_warmup = Nw,
                   num_sample = 0,
                   metric = "dense_e",
                   save_diagnostics = TRUE,
                   inv_metric = inv_metric,
                   init = init_file,
                   stepsize = stepsize,
                   adapt_delta = models[[model_i, "adapt_delta"]],
                   init_buffer = 0,
                   window = Nw + 1,
                   term_buffer = 0)

usamples = getUnconstrainedSamples(fit$diagnostic_files())

plot(1:Nw, usamples[, 3])
Ytrain = head(usamples, -Nw / 2)
Ytest = tail(usamples, Nw / 2)
lpdf = lapply(names(inv_metrics), function(inv_metric_name) {
  print(inv_metric_name)
  inv_metric = inv_metrics[[inv_metric_name]](Ytrain)
  lps = dmvnorm(Ytest, colMeans(Ytrain), inv_metric, log = TRUE)
  neff = monitor(array(lps, dim = c(length(lps), 1, 1)), warmup = 0, print = FALSE)[, "Bulk_ESS"]

  L = t(chol(inv_metric))
  el = eigen(solve(L, t(solve(L, t(cov(Ytest))))), T)
  ew = eigen(solve(L, t(solve(L, t(lw_linear_inv_metric(Ytest))))), T)
  H = t(L) %*% getHessian(stan_fit, tail(usamples, 1)) %*% L
  eh = eigen(H, T)
  out = max(abs(eh$values)) * max(abs(el$values))
  out2 = max(abs(el$values)) / min(abs(el$values))
  out3 = max(abs(eh$values)) * max(abs(ew$values))
  out4 = max(abs(ew$values)) / min(abs(ew$values))

  tibble(name = inv_metric_name,
         c_hybrid = sqrt(out),
         c_hybrid_wolf = sqrt(out3),
         c_matrix = sqrt(out2),
         c_matrix_wolf = sqrt(out4),
         lp = mean(lps),
         sde = sd(lps) / sqrt(neff),
         lp_neff = neff)
}) %>%
  bind_rows %>%
  arrange(lp)

lpdf

name = 
  lpdf %>%
  tail(1) %>%
  pull(name)

init_file = get_init(usamples)
stepsize = get_stepsize(fit)

Ntest = 50
essdf = lapply(names(inv_metrics), function(name) {
  fit = model$sample(data = data,
                     num_chains = 1,
                     save_warmup = 1,
                     num_warmup = Ntest,
                     num_sample = 0,
                     experimental = 0,
                     metric = "dense_e",
                     which_adaptation = 0,
                     save_diagnostics = TRUE,
                     inv_metric = inv_metrics[[name]](Ytrain),
                     init = init_file,
                     stepsize = stepsize,
                     adapt_delta = models[[model_i, "adapt_delta"]],
                     init_buffer = 0,
                     window = Ntest + 1,
                     term_buffer = 0,
                     refresh = 10)
  
  usamples = getUnconstrainedSamples(fit$diagnostic_files())
  monitor(getExtras(fit$diagnostic_files()) %>% pull(lp__) %>% array(dim = c(length(.), 1, 1)), warmup = 0)
  mdf = monitor(array(usamples, dim = c(nrow(usamples), 1, ncol(usamples))), warmup = 0, print = FALSE)
  leapfrogs = getLeapfrogs(fit$diagnostic_files())
  tibble(name = name,
         min_ess = min(mdf[, "Bulk_ESS"]),
         min_ess_p_lf = min_ess / leapfrogs,
         max_rhat = max(mdf[, "Rhat"]))
}) %>% bind_rows()

lpdf %>%
  left_join(essdf)

inv_metric = inv_metrics[[name]](Ytrain)
