########################
# Wrapper for optmatch
########################

library(optmatch)
library(tidyverse)
options("optmatch_max_problem_size" = Inf)

match_fit <- 
  function(
    formula,
    propensity,
    data,
    se_type = "stata",
    clusters = NULL
  ) {
    matches <-
      fullmatch(
        propensity, data = data
      )
    fit <- lm_robust(
      formula,
      fixed_effects = ~ matches,
      data = data,
      se_type = se_type,
      clusters = clusters
    )
    return(fit)
  }