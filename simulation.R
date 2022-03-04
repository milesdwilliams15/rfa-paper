
#######################################
# This file contains the code for the #
# Monte Carlo simulation              #
#######################################


# define parameters -------------------------------------------------------

n = 1000 # n of observations
ate = 5 # size of treatment effect


# define d.g.p. -----------------------------------------------------------

  # the confounder x
x = rnorm(n = n, mean = 50, sd = 10)

# the causal variable z
prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
z = rbinom(n = n, size = 1, prob = prZ)

# the outcome variable y
y = 10 + ate*z + 0.5*x + rnorm(n = n, sd = 10)


# test run with RFA -------------------------------------------------------

# install.packages('devtools')
# devtools::install_github('milesdwilliams15/RFA')
library(RFA)

data = data.frame(y,z,x)

result = rfa(y ~ z + x, data = data)

summary_rfa(result)

# Monte Carlo experiment  -------------------------------------------------

  # Methods compared:
  # 1. OLS w/o control
  # 2. OLS w/ control
  # 3. OLS w/ mean-centered interaction 
  # 4. OLS w/ correct specification
  # 5. P-score matching (logit)
  # 6. P-score matching (logit w/ quad term)
  # 7. GenMatch on x
  # 8. SVA
  # 9. RFA

library(Matching) # for matching
library(e1071)    # for support vector machines
library(randomForest) # for random forests
library(sandwich) # for standard errors and p-values
library(lmtest)   # for coefficient summaries
library(tidyverse)# for grammar
library(gridExtra)

  # iterations
its = 999

  # matrices to fill
results = list()

  # set seed
set.seed(999)

  # the loop
startTime = Sys.time()
for(i in 1:its){
  # simulate d.g.p.
  # the confounder x
  x = rnorm(n = n, mean = 50, sd = 10)
  
  # the causal variable z
  prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
  z = rbinom(n = n, size = 1, prob = prZ)
  
  # the outcome variable y
  y = 1 + ate*z + 0.5*x + x^2 + rnorm(n = n, sd = 10)
  
  # get estimates with methods
  est1 = lm(y ~ z) %>%
    coeftest(., vcov = vcovHC(.))
  est2 = lm(y ~ z + x) %>%
    coeftest(., vcov = vcovHC(.))
  est3 = lm(y ~ z + x + z*I(x - mean(x))) %>%
    coeftest(., vcov = vcovHC(.))
  est4 = lm(y ~ z + x + I(x^2)) %>%
    coeftest(., vcov = vcovHC(.))
  est5 = data.frame(x,y,z) %>%
    mutate(pscore = predict(
      glm(z ~ x, family = binomial),
      type = 'response'
    )) %>%
    Match(Y = .$y, 
          Tr = .$z, 
          M=1, 
          replace=T, 
          X=.$pscore) 
  est6 = data.frame(x,y,z) %>%
    mutate(pscore = predict(
      glm(z ~ x + I(x^2), family = binomial),
      type = 'response'
    )) %>%
    Match(Y = .$y, 
          Tr = .$z, 
          M=1, 
          replace=T, 
          X=.$x) 
  matches = data.frame(x,y,z) %>%
    GenMatch(Tr = .$z,
          M=1,
          X=.$x,
          pop.size=5, max.generations = 10,
          replace = T, print.level = 0)
  matches = matches$matches
  df = list()
  for(j in 1:nrow(matches)){
    df[[j]] = data.frame(y,z)[c(matches[j,1:2]),] %>%
      mutate(match = j, weight = matches[j,3])
  }
  est7 = do.call(rbind,df) %>%
    lm(y ~ z + as.factor(match), data = ., weights = .$weight) %>%
    coeftest(., vcovHC(.))
  est8 = data.frame(x,y,z) %>%
    mutate(
      yres = y - predict(
        svm(y ~ x)
      ),
      zres = z - predict(
        svm(z ~ x)
      )
    ) %>%
    lm(yres ~ zres, data=.) %>%
    coeftest(., vcovHC(.))
  est9 = data.frame(x,y,z) %>%
    mutate(
      yres = y - predict(
        suppressWarnings(randomForest(y ~ x))
      ),
      zres = z - predict(
        suppressWarnings(randomForest(z ~ x))
      )
    ) %>%
    lm(yres ~ zres, data=.) %>%
    coeftest(., vcovHC(.))

  # put results into tidy table
  results[[i]] = tibble(
    term = c(
      'OLS',
      'OLS + control',
      'OLS + interaction',
      'OLS + quadradic',
      'Matching (logit)',
      'Matching (logit + quad)',
      'Matching (genetic)',
      'SVA',
      'RFA'
    ),
    estimate = c(
      est1[2,1],
      est2[2,1],
      est3[2,1],
      est4[2,1],
      est5$est,
      est6$est,
      est7[2,1],
      est8[2,1],
      est9[2,1]
    ),
    std.error = c(
      est1[2,2],
      est2[2,2],
      est3[2,2],
      est4[2,2],
      est5$se.standard,
      est6$se.standard,
      est7[2,2],
      est8[2,2],
      est9[2,2]
    ),
    statistic = c(
      est1[2,3],
      est2[2,3],
      est3[2,3],
      est4[2,3],
      est5$est/est5$se.standard,
      est6$est/est6$se.standard,
      est7[2,3],
      est8[2,3],
      est9[2,3]
    ),
    p.value = c(
      est1[2,4],
      est2[2,4],
      est3[2,4],
      est4[2,4],
      2*pnorm(-abs(est5$est/est5$se.standard)),
      2*pnorm(-abs(est6$est/est6$se.standard)),
      est7[2,4],
      est8[2,4],
      est9[2,4]
    ),
    model = i
  )
  
  # report progress of loop
  cat('Progress: ',round(100*i/its,2),'%\r')
  if(i == its) cat("Done!\n")
  flush.console()
}
endTime = Sys.time()
endTime - startTime

# summarize results -------------------------------------------------------

# Make summary table
smry = do.call(rbind, results) %>%
  group_by(term) %>%
  summarize(
    'MSE of Estimate' = sum(((estimate - ate)/its)^2),
    'Mean Bias' = mean(estimate - ate),
    'Coverage' = mean(ate >= (estimate - 1.96*std.error) &
                        ate <= (estimate + 1.96*std.error))
  ) %>%
  mutate_if(is.numeric, function(x) round(x,2)) %>%
  arrange(`MSE of Estimate`)

# print results
smry

# save in .csv file
write_csv(do.call(rbind, results), 'simulation_output.csv')
write_csv(smry, 'simulation_summary.csv')

# write a stargazer table
stargazer::stargazer(smry, 
                     header = F, 
                     summary = F,
                     title = "Monte Carlo Results"
                     )

do.call(rbind,results) %>%
  ggplot() +
  aes(estimate) + 
  geom_density(fill='grey') +
  geom_vline(xintercept = ate) +
  facet_wrap(~term, scales = 'free_y') +
  labs(
    x = "estimated ATE"
  ) +
  theme_test() + 
  ggsave(
    'plot1.pdf',
    units = 'in',
    height = 4,
    width = 6
  )

do.call(rbind,results) %>%
  filter(term!='OLS',
         term!='OLS + control',
         term!='OLS + interaction') %>%
  ggplot() +
  aes(estimate) + 
  geom_density(fill='grey') +
  geom_vline(xintercept = ate) +
  facet_wrap(~term, scales = 'free_y') +
  labs(
    x = "estimated ATE"
  ) +
  theme_test() + 
  ggsave(
    'plot3.pdf',
    units = 'in',
    height = 3,
    width = 6
  )



# propensity to receive treatment -----------------------------------------

prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
p1 = ggplot() +
  aes(x,prZ) +
  geom_line(size = 1) +
  labs(
    x = expression(x[i]),
    y = expression('Pr'*(z[i]==1*'|'*x[i])),
    caption = "Propensity to receive treatment"
  ) +
  scale_x_continuous(breaks=seq(20,80,by=5)) +
  scale_y_continuous(breaks=c(0,1),
                     limits=c(0,1)) +
  theme_test() 
y = 1 + 0.5*x + x^2
p2 = ggplot() +
  aes(x,y) +
  geom_line(size = 1) +
  labs(
    x = expression(x[i]),
    y = expression('E'*(y[i]*'|'*x[i])),
    caption = "Conditional mean of response"
  ) +
  scale_x_continuous(breaks=seq(20,80,by=5)) +
  theme_test() 
grid.arrange(p1,p2,ncol=2)
ggsave(
  plot = grid.arrange(p1,p2,ncol=2),
  'plot2.pdf',
  units = 'in',
  height = 3, width = 7
)






