library(nlshrink)
library(cmdstanr)
library(tidyverse)
library(rstan)
library(nlshrink)
library(parallel)

source("inv_metrics.R")
set_cmdstan_path("/home/bbales2/warmup/cmdstan/")

models = tibble(name = c("kilpisjarvi", "accel_splines", "accel_gp", "diamonds", "prophet", "radon"),
                adapt_delta = c(0.8, 0.99, 0.99, 0.8, 0.8, 0.8),
                init = c(2.0, 0.25, 0.25, 2.0, 0.5, 2.0))

model_i = 4

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
                   num_chains = 2,
                   save_warmup = 0,
                   experimental = 1,
                   metric = "dense_e",
                   which_adaptation = 0,
                   save_diagnostics = TRUE,
                   init = models[[model_i, "init"]],
                   adapt_delta = models[[model_i, "adapt_delta"]],
                   num_cores = 2)

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

X = getUnconstrainedSamples(fit$diagnostic_files())[1:40,]
D = diag(cov(X))
evecs = eigen(cov(X), T)$vectors
#X = apply(X, 1, function(row) { row / sqrt(D) }) %>% t

N = nrow(X)
P = ncol(X)

#X = sapply(1:N, function(x) rnorm(P)) %>% t



edf = lapply(seq(3, N, by = 1), function(Nm) {
  print(Nm)
  tibble(#lw_evals = eigen(lw_linear_metric(X[1:Nm, ]), T)$values[1:min(P, Nm - 1)],
         #hess_evals = eigen(hessian_metric(X[(Nm - 1):Nm, ]), T)$values[1:min(P, Nm - 1)],
         #lw_corr_evals = eigen(lw_linear_corr_metric(X[1:Nm, ]), T)$values[1:min(P, Nm - 1)],
         samp_evals = eigen(cov(X[1:Nm, ]), T)$values,#[1:min(P, Nm - 1)]
         diag_evals = diag(diag_inv_metric(X[1:Nm, ])),
         dense_evals = eigen(dense_inv_metric(X[1:Nm, ]), T)$values,
         #oracle_evals = diag(t(evecs) %*% cov(X[1:Nm, ]) %*% evecs),
         dense2_evals = eigen(dense2_inv_metric(X[1:Nm, ]), T)$values,
         #eval_bound = sum(diag(diag_inv_metric(X[1:Nm, ]))) / P,
         r = 1:P,#min(P, Nm - 1),
         Nm = Nm)
}) %>% bind_rows

edf %>%
  gather(which, evals, diag_evals, dense_evals, samp_evals, dense2_evals) %>%
  filter(evals > 1e-10) %>%
  ggplot(aes(Nm, evals)) +
  geom_line(aes(color = which, group = interaction(r, which)), size = 0.5, alpha = 0.75, position = position_dodge(width = 1.0)) +
  scale_y_log10()

