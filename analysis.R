
##########################################
# Analysis of Nielsen et al. (2011) data #
##########################################
rm(list = ls())

# load libraries ----------------------------------------------------------

source("source_code/rfa_script.R")
source("source_code/match_fit_script.R")
library(ggridges)
library(stargazer)

# load the dataset --------------------------------------------------------

load("data/civil_war_data.R")
dt <- imp_data         # simpler name
dim(na.omit(raw_data)) - dim(dt) # gained 175 obs.
rm(raw_data, imp_data) # don't need
dt <- dt %>%
  mutate(
    civil_war = ifelse(civil_war < 0.5, 0, 1),
    shock = ifelse(
      aid_change <= quantile(aid_change, 0.15),
      1, 0
    ),
    pos_shock = ifelse(
      aid_change >= quantile(aid_change, 0.85),
      1, 0
    ),
    shock25 = ifelse(
      aid_change <= quantile(aid_change, 0.25),
      1, 0
    ),
    pos_shock75 = ifelse(
      aid_change >= quantile(aid_change, 0.75),
      1, 0
    )
    
  )

# summary stats -----------------------------------------------------------

nrow(dt) # n = 5,946
dt %>% 
  summarize(
    n_countries = length(unique(country))
  )

hist_plot <-
  ggplot(dt) +
  aes(aid_change) +
  geom_density(fill='grey') +
  labs(
    x = "Change in Aid/GDP",
    y = ''
  ) +
  scale_x_continuous(limits = c(-.2,.2)) +
  geom_vline(
    aes(xintercept = quantile(aid_change,.15)),
    lty = 2
  ) +
  theme_ridges() 
hist_plot

dt %>%
  count(civil_war) %>%
  mutate(prop = n/sum(n)) %>%
  stargazer(
    summary = F,
    header = F,
    title = "Instances of Civil War Onset"
  )

dt %>%
  #group_by(neg_shock) %>%
  count(civil_war) %>%
  mutate(prop = n/sum(n)) %>%
  mutate_all(function(x) round(x,2)) %>%
  stargazer(
    summary = F,
    header = F,
    title = "Instances of Civil War Onset"
  )
# library(foreach)
# cis = foreach(i = 1:1000, .combine = 'rbind') %dopar% {
#   require(dplyr)
#   dt %>%
#     sample_n(size = nrow(dt), replace = T) %>%
#     group_by(neg_shock) %>%
#     count(civil_war) %>%
#     mutate(prop = n/sum(n)) %>%
#     filter(civil_war == 1) %>%
#     .$prop
# } %>%
#   apply(.,2,function(x) quantile(x, c(.025,.975))) %>%
#   t(.)
# 
# dt %>%
#   group_by(neg_shock) %>%
#   count(civil_war) %>%
#   mutate(prop = n/sum(n)) %>%
#   filter(civil_war == 1) %>%
#   ungroup() %>%
#   mutate(
#     lo.prop = cis[,1],
#     hi.prop = cis[,2]
#   ) %>%
#   ggplot() + 
#   aes(neg_shock,prop,
#       ymin = lo.prop,
#       ymax = hi.prop) +
#   geom_col(width = .5) +
#   geom_errorbar(width = .15) +
#   scale_x_discrete(
#     labels = c("No Shock","Shock")
#   ) +
#   labs(
#     x = "",
#     y = "Proportion"
#   ) +
#   theme_ridges() +
#   theme(
#     text=element_text(family='Palatino Linotype')
#   ) +
#   coord_flip() +
#   ggsave(
#     "diffprop.png",
#     height = 3,
#     width = 6
#   )


# analysis ----------------------------------------------------------------


 # Main results
rfa_fit_neg <- rfa(
  civil_war ~ shock,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)
rfa_fit_neg25 <- rfa(
  civil_war ~ shock25,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)
rfa_fit_pos <- rfa(
  civil_war ~ pos_shock,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)
rfa_fit_pos75 <- rfa(
  civil_war ~ pos_shock75,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)

 # Compare results
texreg::texreg(
  list(rfa_fit_neg$fit,
       rfa_fit_pos$fit,
       rfa_fit_neg25$fit,
       rfa_fit_pos75$fit),
  include.ci = F,
  digits = 3,
  custom.model.names = 
    c("Negative", "Positive",
      "Negative (25th)", "Positive (75th)"),
  custom.coef.names = 
    c("Control", "Shock"),
  include.rsquared = F,
  include.adjrs = F,
  include.rmse = F,
  caption = "Random Forest Adjusted Estimates",
  caption.above = T
)

# Run loop with different cutpoints for aid shocks.

# RFA analysis first:

pcts <- seq(.2,.8,by=.05)
pcts <- pcts[which(pcts!=.5)]
out <- list()
for(i in 1:length(pcts)) {
  if(i == 1) cat("Starting",
                 rep(".",len=10),
                 "|\n")
  if(pcts[i]<.5){
    new_dt <- dt %>%
      mutate(
        shock = 
          ifelse(
            aid_change <= quantile(aid_change,pcts[i]),
            1, 0
          )
      )
  } else {
    new_dt <- dt %>%
      mutate(
        shock = 
          ifelse(
            aid_change >= quantile(aid_change,pcts[i]),
            1, 0
          )
      )
  }
  
  # RFA
  rfa_fit <- suppressMessages(rfa(
    civil_war ~ shock,
    ~ pos_shock + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
      x20 + year,
    new_dt,
    clusters = new_dt$country
  )$fit %>% tidy())
  # ipw_fit <- rfa_rose(
  #   civil_war ~ shock,
  #   ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
  #     x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
  #     x20 + as.factor(year),
  #   new_dt,
  #   clusters = new_dt$country
  # ) %>% tidy()
  # ols_fit <- suppressMessages(lm_robust(
  #   civil_war ~ shock + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
  #     x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
  #     x20 + as.factor(year),
  #   new_dt,
  #   clusters = new_dt$country,
  #   se_type = "stata"
  # ) %>% tidy())
  # # lin_fit <- suppressMessages(lm_lin(
  # #   civil_war ~ shock,
  # #   ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
  # #     x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
  # #     x20 + as.factor(year),
  # #   new_dt,
  # #   clusters = new_dt$country,
  # #   se_type = "stata"
  # # ) %>% tidy())
  # mch_fit <- match_fit(
  #   civil_war ~ shock,
  #   shock ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
  #     x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 +
  #     x20 + year,
  #   new_dt,
  #   clusters = new_dt$country
  # ) %>% tidy()
  out[[i]] <- bind_rows(
    rfa_fit
    # ipw_fit,
    # ols_fit,
    # # lin_fit,
    # mch_fit
  ) %>%
    filter(
      term %in% c("xres", "shock")
    ) %>% 
    mutate(
      # model = c("RFA", "MRA", "Matching"),
      cutoff = pcts[i]
    ) 
  cat("Progress",
      rep(".",len=round(10 * i / length(pcts))),
      rep(" ",len=10-round(10*i / length(pcts))),
      "|\r\r\r\r\r\r\r\r\r\r\r\r")
  if(i == length(pcts)) cat("\nDone!")
} 

results <-
  bind_rows(out)

civ_plot <- results %>%
  ggplot() +
  aes(
    x = cutoff,
    y = estimate,
    ymin = conf.low,
    ymax = conf.high
  ) +
  geom_line() +
  geom_ribbon(
    alpha = 0.5
  ) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(
    breaks = pcts,
  ) +
  labs(
    x = "Percentile Cutoff",
    y = "Estimated Effect with 95% CI",
    color = NULL
  ) +
  annotate(
    "text",
    x = c(.35,.65),
    y = c(-.025,.025),
    label = c(
      "''%<-%'Neg. Shock'",
      "'Pos. Shock'%->%''"
    ),
    parse = T
  ) +
  theme_ridges(center_axis_labels = T) + 
  theme(
    legend.position = 'top',
    strip.text = element_blank()
  ) 
save(civ_plot, file = "figures/aid_shock_effect.R")

#### With a continuous measure of aid shock
stand <- function(x) (x - mean(x)) / sd(x)
dt <- dt %>% 
  mutate(
    aid_quants = gtools::quantcut(
      aid_change, q = 3
    ) %>% as.numeric() - 2,
    aid_rank = stand(rank(aid_change)),
    aid_shock = ifelse(
      aid_change <= quantile(aid_change, 0.25),
      1, 0
    )
  )

# RFA
rfa_quant <- rfa(
  civil_war ~ aid_quants,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)
rfa_change <- rfa(
  civil_war ~ aid_change,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)
rfa_rank <- rfa(
  civil_war ~ aid_rank,
  ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +
    x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + 
    x20 + year,
  dt,
  clusters = dt$country,
  importance = "impurity"
)
texreg::texreg(
  list(
    rfa_change$fit,
    rfa_rank$fit,
    rfa_quant$fit
  ),
  custom.model.names = 
    c("Aid/GDP", "Rank", "Percentile"),
  custom.coef.names = 
    c("Constant","Change"),
  include.ci = F,
  include.rsquared = F,
  include.adjrs = F,
  include.rmse = F,
  digits = 3,
  caption = "Random Forest Adjusted Estimates (Continuous Predictor)",
  caption.above = T
)


# Permutation importance of confounders -----------------------------------

impy <- importance(rfa_fit_neg25$yrf)
impx <- importance(rfa_fit_neg25$xrf)
names(impy) <- names(impx) <- c(covnames, "Year")
impdf <- bind_cols(
  var = c(covnames, "Year"),
  "Civil War" = impy,
  "Aid Shock" = impx
) %>%
  gather(
    key = "value",
    value = "Importance",
    -var
  )
rho <- cor(impy, impx)
ggplot(impdf) +
  aes(
    Importance,
    reorder(var, Importance),
    color = value
  ) +
  geom_vline(
    xintercept = 0,
    alpha = 0
  ) +
  geom_point() +
  labs(
    y = NULL,
    color = NULL,
    caption = 
      substitute(
        paste(rho, " = ", r), list(r = round(rho, 2))
      )
  ) +
  theme_ridges(
    font_size = 10,
    center_axis_labels = T
  ) +
  theme(
    legend.position = "top"
  ) +
  ggsave(
    "figures/varimport.png",
    height = 6,
    width = 6
  )


# Nominal vs. effective sample --------------------------------------------

w <-
  rfa_fit_neg25$data$xres^2
wmean <- 
  function(x) {
    wmean <- sum(w * x) / sum(w)
    return(wmean)
  }
wse <- function(x) {
  wse <- diagis::weighted_se(
    x, w = w
  )
  return(wse)
}

samp_means <- 
  dt %>%
  mutate_if(
    is.numeric,
    function(x) (x - mean(x)) / sd(x)
  ) %>%
  summarize_if(
    is.numeric, 
    list(wmean, wse)
  ) %>%
  select(
    contains("x")
  ) %>%
  gather() %>%
  mutate(
    vars = c(covnames, covnames),
    key = rep(c("est", "se"), each = 20)
  ) %>%
  spread(
    key = key,
    value = value
  )

ggplot(samp_means) +
  aes(
    est,
    reorder(vars, est^2),
    xmin = est - 1.96 * se,
    xmax = est + 1.96 * se
  ) +
  geom_vline(
    xintercept = 0,
    alpha = 0
  ) +
  geom_point() +
  geom_errorbarh(
    height = 0
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  labs(
    x = "Standard Deviation Units",
    y = NULL,
    caption = "0 = Nominal Sample Mean"
  ) +
  theme_ridges(
    center_axis_labels = T,
    font_size = 10
  ) +
  ggsave(
    "figures/effsampdiff.png",
    height = 6,
    width = 6
  )