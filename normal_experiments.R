library(tidyverse)
library(cmdstanr)
library(rstan)

set_cmdstan_path("/home/bbales2/tmp/cmdstan")

model = cmdstan_model("normal.stan", quiet = FALSE)

data = list()
for(N in seq(40, 400, 40)) {
  for(metric in c("dense_e", "diag_e")) {
    out = model$sample(data = list(N = N),
                       adapt_engaged = FALSE,
                       num_warmup = 0,
                       num_samples = 2000,
                       init = 0,
                       seed = 2020,
                       metric = metric,
                       num_chains = 4)
    
    fit = read_stan_csv(out$save_output_files())

    leapfrogs = get_num_leapfrog_per_iteration(fit) %>% sum
    time = sum(get_elapsed_time(fit)[, "sample"])
    print(paste0("N: ", N, ", metric: ", metric, ", time :", time))
    data[[length(data) + 1]] = tibble(time = time / leapfrogs,
                                      N = N,
                                      metric = metric)
  }
}

df = data %>%
  bind_rows() %>%
  spread(metric, time) %>%
  mutate(ratio = dense_e / diag_e)

df %>%
  ggplot(aes(N, ratio)) +
  geom_point(aes(color = metric)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
