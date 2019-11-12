library(tidyverse)
library(cmdstanr)
library(rstan)

cmdstan_path = "/home/bbales2/cmdstan-warmup/"

set_cmdstan_path(cmdstan_path)

models = tibble(name = c("kilpisjarvi", "accel_splines", "accel_gp", "diamonds", "prophet", "radon"),
                adapt_delta = c(0.8, 0.99, 0.99, 0.8, 0.8, 0.8),
                init = c(2.0, 0.25, 0.25, 2.0, 2.0, 2.0))

num_repetitions = 4
num_chains = 4
metric = "dense_e"

initdf = lapply(1:nrow(models), function(i) {
  model_name = models[[i, "name"]]
  
  model = cmdstan_model(paste0(cmdstan_path, "examples/", model_name, "/", model_name, ".stan"), quiet = FALSE)
  
  data_env = new.env(parent = baseenv())
  source(paste0(cmdstan_path, "examples/", model_name, "/", model_name, ".dat"), local = data_env)
  
  # Find typical set
  fit <- model$sample(data = as.list.environment(data_env),
                      num_chains = 1,
                      experimental = 1,
                      num_samples = 0,
                      num_warmup = 150,
                      which_adaptation = 0,
                      adapt_delta = models[[i, "adapt_delta"]],
                      init = models[[i, "init"]],
                      metric = "dense_e",
                      save_warmup = TRUE,
                      save_diagnostics = TRUE,
                      init_buffer = 150,
                      window = 0,
                      term_buffer = 0)
  
  diagdf = read_csv(fit$diagnostic_files(), comment = "#") %>%
    mutate(draw = row_number())
  
  nparams = ncol(diagdf %>% select(starts_with("p_")))
  
 p = diagdf %>%
    mutate(momentum_squared = rowSums((diagdf %>%
                                          select(starts_with("p_")) %>%
                                          as.matrix())^2),
           sum_gradients = rowSums(diagdf %>%
                                     select(starts_with("g_")) %>%
                                     as.matrix())) %>%
    mutate(model = model_name) %>%
    select(draw, model, momentum_squared, sum_gradients, lp__, divergent__) %>%
    gather(variable, value, -model, -draw, -divergent__) %>%
    group_by(variable, model) %>%
    mutate(group_mean = mean(tail(value, 50)),
           group_sd = sd(tail(value, 50))) %>%
    mutate(values_off_plot = abs((value - group_mean) / group_sd) > 20.0,
           value = ifelse(values_off_plot, sign(value) * 20 * group_sd + group_mean, value)) %>%
    ggplot(aes(draw, value)) +
    geom_point(aes(color = as.factor(values_off_plot)), alpha = 0.75) +
    facet_grid(variable ~ ., scales = "free_y") +
   ggtitle(paste0("Model: ", model_name, ", # params: ", nparams, ", far away values truncated"))
  
  ggsave(paste0("typical_set/", model_name, ".png"), plot = p, scale = 1, width = 9, height = 5, dpi = 300)
}) %>% bind_rows

initdf %>%
  filter(model == "diamonds") %>%
  select(r, model, scaled_p_mag, sum_gradients, lp__, divergent__) %>%
  gather(variable, value, -model, -r, -divergent__) %>%
  group_by(variable, model) %>%
  mutate(value = (value - mean(tail(value, 50))) / sd(tail(value, 50))) %>%
  filter(abs(value) < 10000.0) %>%
  ggplot(aes(r, value)) +
  geom_point(aes(color = model, shape = as.factor(divergent__)), alpha = 0.75) +
  facet_grid(variable ~ ., scales = "free_y")
