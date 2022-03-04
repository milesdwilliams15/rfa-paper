##############################################
# Script for RFA function plus other helpers
##############################################


# Required packages -------------------------------------------------------

library(tidyverse) # For grammar
library(fabricatr) # For simulating data
library(ranger)    # For fast rf implementation
library(estimatr)  # For linear model fitting and inference


# Function to simulate data -----------------------------------------------

sim_data <- function(N, ate = 0, nl = TRUE) {
  if(nl == TRUE) {
    fabricate(
      N = N,
      x = rnorm(N, 50, sd = 10),
      prZ = 1 / 
        (1 + exp(2 - (0.05 * (1.15 * x - mean(x))^2))),
      z = rbinom(N, 1, prZ),
      y = 1 + ate * z + 0.5 * x + x^2 + rnorm(N, sd = 10)
    )
  } else {
    fabricate(
      N = N,
      x = rnorm(N, 50, sd = 10),
      prZ = 1 / 
        (1 + exp(2 - (0.05 * (1.15 * x)))),
      z = rbinom(N, 1, prZ),
      y = 1 + ate * z + 0.5 * x + rnorm(N, sd = 10)
    )
  }
}


# rfa ---------------------------------------------------------------------

rfa <-
  function(
    formula,
    covariates,
    data,
    se_type = "stata",
    clusters = NULL,
    ...
  ) {
    
    # Get design matrix of covariates
    covmat <-
      model.matrix(
        covariates, data = data
      )[, -1]
    
    # Vectors of treatment and response 
    varmat <-
      cbind(
        model.frame(
          formula, data = data
        )[, 1],
        model.matrix(
          formula, data = data
        )[, -1]
      )
    
    # Residulalize the response
    yrf <-
      ranger(
        y = varmat[, 1],
        x = cbind(covmat),
        # probability = ifelse(
        #   length(unique(varmat[, 1])) == 2,
        #   TRUE, FALSE
        # ),
        ...
      )
    yhat <- yrf$predictions
    if(!is.null(dim(yhat))) yhat <- yhat[, 2]
    yres <- varmat[, 1] - yhat
    
    # Residualize the treatment
    xrf <-
      ranger(
        y = varmat[, 2],
        x = cbind(covmat),
        # probability = ifelse(
        #   length(unique(varmat[, 2])) == 2,
        #   TRUE, FALSE
        # ),
        ...
      )
    xhat <- xrf$predictions
    if(!is.null(dim(xhat))) xhat <- xhat[, 2]
    xres <- varmat[, 2] - xhat
    
    # Estimate
    data <-
      data %>%
      mutate(yres = yres, xres = xres)
    fit <-
      lm_robust(
        yres ~ xres,
        # update(covariates, yres ~ xres + .),
        data = data,
        se_type = se_type,
        clusters = clusters
      )
    
    # Return fitted model
    lst <- list(
      fit = fit,
      yrf = yrf,
      xrf = xrf,
      data = data
    )
    return(lst)
  }


# rfa_rose ----------------------------------------------------------------

rfa_rose <-
  function(
    formula,
    covariates,
    data,
    se_type = "stata",
    clusters = NULL,
    ...
  ) {
    
    # Get design matrix of covariates
    covmat <-
      model.matrix(
        covariates, data = data
      )[, -1]
    
    # Vectors of treatment and response 
    varmat <-
      cbind(
        model.frame(
          formula, data = data
        )[, 1],
        model.matrix(
          formula, data = data
        )[, -1]
      )
    
    # Residulalize the response
    yhat <-
      ranger(
        y = varmat[, 1],
        x = cbind(covmat),
        # probability = ifelse(
        #   length(unique(varmat[, 1])) == 2,
        #   TRUE, FALSE
        # ),
        ...
      )$predictions
    if(!is.null(dim(yhat))) yhat <- yhat[, 2]
    yres <- varmat[, 1] - yhat
    
    # Get inverse probability weights for treatment
    xhat <-
      ranger(
        y = varmat[, 2],
        x = cbind(covmat),
        # probability = ifelse(
        #   length(unique(varmat[, 2])) == 2,
        #   TRUE, FALSE
        # ),
        ...
      )$predictions 
    if(!is.null(dim(xhat))) xhat <- xhat[, 2]
    
    # Get inverse probability weights:
    if(length(unique)==2) {
      # For a binary predictor:
      ipw <- varmat[, 2] * (1 / xhat) +
        (1 - varmat[, 2]) * (1 / (1 - xhat))
    } else {
      # For a continuous predictor:
      # see from Naimi et al. (2014) on Andrew Heiss' blog
      ipw <- 
        dnorm(
          varmat[, 2],
          mean(varmat[, 2]),
          sd(varmat[, 2])
        ) /
        dnorm(
          varmat[, 2],
          xhat,
          sd(varmat[, 2] - xhat)
        )
    }
    
    # Estimate
    data <-
      data %>%
      mutate(yres = yres, ipw = ipw)
    fit <-
      lm_robust(
        update(formula, yres ~ .),
        data = data,
        weights = ipw,
        se_type = se_type,
        clusters = clusters
      )
    
    # Return fitted model
    return(fit)
  }

# # Test
# rfa(
#   voted ~ got,
#   covariates = ~ female + income + education + partisan,
#   data = sim_data()
# ) %>% tidy()

# rfa_rose(
#   voted ~ got,
#   covariates = ~ female + income + education + partisan,
#   data = df
# )
