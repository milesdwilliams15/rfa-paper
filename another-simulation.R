################################
# Use Sharp Null to Refine RFA #
################################



# libraries ---------------------------------------------------------------

library(tidyverse)
library(randomForest)

# rfa_model ---------------------------------------------------------------

rfa_model <- function(
  form,
  data=NULL,
  ntree=c(500,500),
  replace=c(T,T),
  sampsize=c(0.632,0.632),
  boot=10000
) {
  
  # get data
  data <- model.frame(form,data=data)
  
  # rename response and treatment
  colnames(data) <- c('y','z',paste0('x',1:(ncol(data)-2)))
  
  # residualize response and treatment as a function
  # of covariates
  yhat <- suppressWarnings(
    predict(
      randomForest(
        y ~ ., data = data[,-2],
        ntree = ntree[1],
        replace = replace[1],
        sampsize = ceiling(sampsize[1]*nrow(data))
      )
    )
  )
  zhat <- suppressWarnings(
    predict(
      randomForest(
        z ~ ., data = data[,-1],
        ntree = ntree[2],
        replace = replace[2],
        sampsize = ceiling(sampsize[2]*nrow(data))
      )
    )
  )
  
  yres <- data[,1] - yhat
  zres <- data[,2] - zhat
  
  # estimate ATE
  ATE <- coef(lm(yres~zres))[2]
  
  # bootstrap to get std. error
  its <- boot
  boots <- 0
  for(i in 1:its){
    cat('Bootstrapping..........',i,'of',its,'\r\r\r\r\r')
    obs <- sample(1:nrow(data),size=nrow(data),
                  replace=T)
    boots[i] <- coef(lm(yres[obs]~zres[obs]))[2]
    if(i==its) cat('\nDONE!')
  }
  
  # return ATE and std.error + replicates
  return(
    list(
      results = data.frame(
        estimate = ATE,
        std.error = sd(boots)
      ),
      boots = boots,
      model = data
    )
  )
}

tidy_rfa <- function(x,model) as_tibble(
    x$results
  ) %>%
    mutate(
      term = 'treatment',
      statistic = estimate/std.error,
      p.value = 2*pnorm(abs(statistic)),
      conf.low = estimate - 1.96*std.error,
      conf.high = estimate + 1.96*std.error,
      model = model
    ) %>%
    select(
      term, everything()
    )


# simulate data -----------------------------------------------------------

N <- 10000
tibble(
  x1 = rnorm(N),
  x2 = rbinom(N,1,prob=0.2),
  y0.p = pnorm(
    -1 - x1 + .1*x1^2 + 0.01*x2
  ),
  y1.p = pnorm(
    -1 - x1 + .1*x1^2 + 0.01*x2 - 0.05
  ),
  tr.p = pnorm(
    0.1 + x1 - .2*x1^2 + 0.01*x2
  ),
  tr.cond = rbinom(N,1,tr.p),
  y.cond = rbinom(N,1,y1.p*tr.cond+y0.p*(1-tr.cond)),
  tr.rand = rbinom(N,1,.5),
  y.rand = rbinom(N,1,y1.p*tr.rand+y0.p*(1-tr.rand))
) -> data


# comparison to Niave, Additive Controls, Lin -----------------------------

library(estimatr)

niave_ols <- lm_robust(
  y.rand ~ tr.rand, data = data
)
addit_ols <- lm_robust(
  y.rand ~ tr.rand + x1 + x2, data = data
)
lin_ols <- lm_lin(
  y.rand ~ tr.rand, covariates = ~ x1 + x2,
  data = data
)
rfa_ols <- rfa_model(
  y.rand ~ tr.rand + x1 + x2, data = data,
  boot = 100
)

tidy_models <- function(x,model) tidy(x) %>%
  select(-df,-outcome) %>%
  .[2,] %>%
  mutate(model = model)

rbind(
  tidy_models(niave_ols,'Niave'),
  tidy_models(addit_ols,'Additive'),
  tidy_models(lin_ols,'Lin'),
  tidy_rfa(rfa_ols,'RFA')
) -> smry

ggplot(smry) +
  aes(
    reorder(model,-abs(estimate-0.05)),
    estimate,
    ymin=conf.low,
    ymax=conf.high
  ) +
  geom_point() +
  geom_errorbar(width=0) +
  geom_hline(yintercept = -0.05) +
  labs(
    x = '',
    y = 'ATE (with 95% CIs)'
  )

