########################
# Monte Carlo analysis 
########################
rm(list = ls())


# Attach source code ------------------------------------------------------

source("source_code/rfa_script.R")
source("source_code/match_fit_script.R")


# Run simulation ----------------------------------------------------------

 # Simulation set-up:
sim_output_nl <- sim_output_l <- list()
nsims <- 100
nsamp <- c(500, 1000, 2000)
grid_search <- expand.grid(
  sim = 1:nsims,
  nsamp = nsamp
)
set.seed(999)

 # Nonlinear confounding:
for(i in 1:nrow(grid_search)) {
  if(i == 1) cat("Working on nonlinear confounding ........\n")
  true_ate <- 5
  df <- sim_data(grid_search$nsamp[i], ate = true_ate)
  naive_ate <- with(df, mean(y[z==1]) - mean(y[z==0]))
  rfa_fit <- rfa(
    y ~ z,
    covariates = ~ x,
    data = df
  )$fit %>% tidy() 
  # ros_fit <- rfa_rose(
  #   y ~ z,
  #   covariates = ~ x,
  #   data = df
  # ) %>% tidy()
  ols_fit <- lm_robust(
    y ~ z + x,
    data = df,
    se_type = "stata"
  ) %>% tidy()
  # lin_fit <- lm_lin(
  #   y ~ z,
  #   covariates = ~ x,
  #   data = df,
  #   se_type = "stata"
  # ) %>% tidy()
  # mch_fit <- match_fit(
  #   y ~ z,
  #   propensity = z ~ x,
  #   data = df
  # ) %>% tidy()
  sim_output_nl[[i]] <- 
    bind_rows(
      rfa_fit,
      # ros_fit,
      ols_fit
      # lin_fit,
      # mch_fit
    ) %>%
    filter(
      term %in% c("xres", "z")
    ) %>%
    mutate(
      estimator = c("RFA", "MRA"),
      true_ate = true_ate,
      naive_ate = naive_ate,
      N = grid_search$nsamp[i]
    )
  cat("% Progress:", 
      round(100 * (i / nrow(grid_search)), 2),
      "..........\r\r\r\r\r\r")
  if(i == nrow(grid_search)) cat("\nDone!")
  rm(df, rfa_fit, ols_fit)
}

 # Linear confounding:
for(i in 1:nrow(grid_search)) {
  if(i == 1) cat("Working on linear confounding ........\n")
  true_ate <- 5
  df <- sim_data(grid_search$nsamp[i], ate = true_ate,
                 nl = FALSE)
  naive_ate <- with(df, mean(y[z==1]) - mean(y[z==0]))
  rfa_fit <- rfa(
    y ~ z,
    covariates = ~ x,
    data = df
  )$fit %>% tidy() 
  # ros_fit <- rfa_rose(
  #   y ~ z,
  #   covariates = ~ x,
  #   data = df
  # ) %>% tidy()
  ols_fit <- lm_robust(
    y ~ z + x,
    data = df,
    se_type = "stata"
  ) %>% tidy()
  # lin_fit <- lm_lin(
  #   y ~ z,
  #   covariates = ~ x,
  #   data = df,
  #   se_type = "stata"
  # ) %>% tidy()
  # mch_fit <- match_fit(
  #   y ~ z,
  #   propensity = z ~ x,
  #   data = df
  # ) %>% tidy()
  sim_output_l[[i]] <- 
    bind_rows(
      rfa_fit,
      # ros_fit,
      ols_fit
      # lin_fit,
      # mch_fit
    ) %>%
    filter(
      term %in% c("xres", "z")
    ) %>%
    mutate(
      estimator = c("RFA","MRA"),
      true_ate = true_ate,
      naive_ate = naive_ate,
      N = grid_search$nsamp[i]
    )
  cat("% Progress:", 
      round(100 * (i / nrow(grid_search)), 2),
      "..........\r\r\r\r\r\r")
  if(i == nrow(grid_search)) cat("\nDone!")
  rm(df, rfa_fit, ols_fit)
}

sim_output <- 
  bind_rows(
    bind_rows(
      sim_output_nl
    ) %>% mutate(confounding = "Nonlinear"),
    bind_rows(
      sim_output_l
    ) %>% mutate(confounding = "Linear")
  )

smry <- sim_output %>%
  group_by(estimator, confounding, N) %>%
  summarize(
    "% Bias" = 100 * mean((estimate - true_ate)/true_ate),
    RMSE = sqrt(mean((estimate - true_ate)^2)),
    "% Coverage" = 100 * mean(conf.low < true_ate & conf.high > true_ate),
    #"% Power" = 100 *mean(p.value <= 0.05)
    "Mean S.E." = mean(std.error)
  ) %>%
  gather(
    key = "metric",
    value = "value",
    -estimator, -confounding, -N
  ) 

g1 <- ggplot(smry) +
  aes(
    x = N,
    y = value,
    color = estimator
  ) +
  geom_line(
    size = 1
  ) +
  scale_x_log10(
    breaks = c(500, 1000, 2000),
    labels = scales::comma
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  facet_wrap(
    metric~confounding, scales = "free",
    ncol = 2
  ) +
  labs(
    y = NULL,
    color = NULL
  ) +
  ggridges::theme_ridges(
    center_axis_labels = T,
    font_size = 10
  ) +
  theme(
    legend.position = "top"
  )
save(g1, file = "figures/monte_carlo_results.R")
