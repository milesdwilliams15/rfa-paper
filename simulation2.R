
#################################
# Nominal vs. Effective Samples #
#################################



# libraries ---------------------------------------------------------------

library(RFA)
library(randomForest)
library(MASS)
library(RItools)
library(optmatch)
library(tidyverse)
library(ggridges)
library(foreach)
library(doParallel)
library(gridExtra)
library(lmtest)
library(sandwich)
extrafont::loadfonts(device='win',quiet=T)

# functions ---------------------------------------------------------------

sample.mean = function(x,w){
  nom_mean = mean(x)
  eff_mean = sum(w*x)/sum(w)
  times(999) %do% {
    obs = sample(1:length(x),replace=T,size=length(x))
    mean(x[obs])
  } -> boot_nom
  times(999) %do% {
    obs = sample(1:length(x),replace=T,size=length(x))
    sum(x[obs]*w[obs])/sum(w[obs])
  } -> boot_eff
  means = c(nom_mean,eff_mean)
  names(means) = c('nominal','effective')
  boot_means = tibble(
    nominal = boot_nom,
    effective = boot_eff
  )
  return(
    list(
      sample_means = means,
      boot = boot_means
    )
  )
}

sample.test = function(x,w){
  w = w/sum(w)
  diff = sum(w*x)/sum(w) - mean(x) 
  boot = times(999) %do% {
    obs = sample(1:length(x),size=length(x),replace=T)
    sum(w[obs]*x[obs])/sum(w[obs]) - mean(x[obs])
  }
  null = boot - diff
  p.value = 2*min(mean(null>=diff),mean(null<=diff))
  return(list(diff = diff,boot = boot, p.value = p.value))
}

overlap = function(x,w){
  require(foreach)
  d1 = density(x,from=min(c(x,w*x)),
               to=max(c(x,w*x)))
  d2 = density(x,from=min(c(x,w*x)),
               to=max(c(x,w*x)),
               weights=w/sum(w))
  overlap = sum(
    apply(rbind(d1$y/sum(d1$y),d2$y/sum(d2$y)),2,min)
  )
  return(overlap)
}


# simulate data -----------------------------------------------------------

set.seed(111)

n = 1000 # n of observations
ate = 5 # size of treatment effect

# the confounder x
x = rnorm(n = n, mean = 50, sd = 10)

# the causal variable z
prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
z = rbinom(n = n, size = 1, prob = prZ)

# the outcome variable y
y = 1 + ate*z + 0.5*x + x^2 + rnorm(n = n, sd = 10)

data = data.frame(y,z,x)

ggplot() +
  aes(
    x,prZ
  ) +
  geom_line(size=1) +
  scale_y_continuous(breaks=seq(0,1,by=.2),
                     limits=c(0,1)) +
  scale_x_continuous(breaks=seq(20,80,by=10)) +
  labs(
    x = expression(x[i]),
    y = 'Probability of Treatment'
  ) +
  theme_ridges() + 
  theme(
    text = element_text(family='Palatino Linotype')
  )-> p1
ggplot() +
  aes(
    x,
    1 + .5*x + x^2
  ) +
  geom_line(size=1) +
  scale_y_continuous(
    breaks=seq(0,6600,by=500),
    labels=scales::comma_format()
  ) +
  scale_x_continuous(breaks=seq(20,80,by=10)) +
  labs(
    x = expression(x[i]),
    y='Expected Value of Response'
  ) +
  theme_ridges() +
  theme(
    text = element_text(family='Palatino Linotype')
  )-> p2
grid.arrange(p1,p2,ncol=2)
ggsave(plot=grid.arrange(p1,p2,ncol=2),
       'treatment.png',
       height=4,
       width=8)

# get weights -------------------------------------------------------------

# matches
matches = fullmatch(z~x,data=data)

# weights
w_ols = resid(lm(z~x,data))^2
w_mtc = resid(lm(z~matches))^2
w_rfa = (
  z - predict(randomForest(z~x))
)^2


# sample differences and overlap ------------------------------------------


test_ols = sample.test(x=x,w=w_ols)
test_mtc = sample.test(x=x,w=w_mtc)
test_rfa = sample.test(x=x,w=w_rfa)

mean_ols = sample.mean(x=x,w=w_ols)
mean_mtc = sample.mean(x=x,w=w_mtc)
mean_rfa = sample.mean(x=x,w=w_rfa)

diff_smry = tibble(
  Estimator = c('OLS','Matching','RFA'),
  'Nominal Mean' = c(round(mean(x),2),NA,NA),
  'Effective Mean' = c(
    sum(w_ols*x)/sum(w_ols),
    sum(w_mtc*x)/sum(w_mtc),
    sum(w_rfa*x)/sum(w_rfa)
  ),
  Difference = c(test_ols$diff,
                 test_mtc$diff,test_rfa$diff),
  'p-value' = c(test_ols$p.value,test_mtc$p.value,
              test_rfa$p.value)
) %>% mutate_if(is.numeric,function(x)round(x,2))
diff_smry
stargazer::stargazer(
  diff_smry,
  summary=F,
  title='Nominal vs. Effective Sample Means',
  header=F,
  rownames=F
)

tibble(
  method = c(
    rep('Regression',2),
    rep('Matching',2),
    rep('RFA',2)
  ),
  sample = c(
    rep(c('nominal','effective'),len=6)
  ),
  mean = c(
    mean_ols$sample_means[1],
    mean_ols$sample_means[2],
    mean_mtc$sample_means[1],
    mean_mtc$sample_means[2],
    mean_rfa$sample_means[1],
    mean_rfa$sample_means[2]
  ),
  lo = c(
    quantile(mean_ols$boot$nominal,.025),
    quantile(mean_ols$boot$effective,.025),
    quantile(mean_mtc$boot$nominal,.025),
    quantile(mean_mtc$boot$effective,.025),
    quantile(mean_rfa$boot$nominal,.025),
    quantile(mean_rfa$boot$effective,.025)
  ),
  hi = c(
    quantile(mean_ols$boot$nominal,.975),
    quantile(mean_ols$boot$effective,.975),
    quantile(mean_mtc$boot$nominal,.975),
    quantile(mean_mtc$boot$effective,.975),
    quantile(mean_rfa$boot$nominal,.975),
    quantile(mean_rfa$boot$effective,.975)
  )
) %>%
  ggplot() +
  aes(
    x = mean,
    y = reorder(method,-mean),
    xmin = lo,
    xmax = hi,
    color = sample
  ) +
  geom_point() +
  geom_errorbarh(height = 0.15) +
  geom_vline(xintercept=mean(x),lty=2) +
  labs(
    x='Effective vs. Nominal Sample Mean',
    y='',
    title='Confounding Variable Difference between Samples'
  ) +
  theme_ridges(
    font_family='Palatino Linotype'
  ) + 
  theme(
    plot.title.position = 'plot',
    legend.title = element_blank(),
    legend.position = 'bottom'
  )

ovl_ols = overlap(x=x,w=w_ols)
ovl_mtc = overlap(x=x,w=w_mtc)
ovl_rfa = overlap(x=x,w=w_rfa)

ovl_smry = tibble(
  estimator = c('OLS','Matching','RFA'),
  '% overlap' = c(ovl_ols,ovl_mtc,ovl_rfa)*100
)
ovl_smry

p1 = ggplot() +
  geom_density(aes(x),color='red',fill='red',alpha=.5) +
  geom_density(aes(sjstats::weight(x,w_ols)),
               color='blue',fill='blue',alpha=.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        text=element_text(family='Palatino Linotype')) +
  labs(
    x='Covariate',
    y='Density',
    title='OLS'
  ) +
  annotate(
    'text',
    x=50,
    y = .06,
    label='Nominal',
    color='red',
    family='Palatino Linotype'
  ) +
  annotate(
    'text',
    x=20,
    y=.055,
    label='Effective',
    color='blue',
    family='Palatino Linotype'
  )

p2 = ggplot() +
  geom_density(aes(x),color='red',fill='red',alpha=.5) +
  geom_density(aes(sjstats::weight(x,w_mtc)),
               color='blue',fill='blue',alpha=.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        text=element_text(family='Palatino Linotype')) +
  labs(
    x='Covariate',
    y='Density',
    title='Matching'
  ) +
  annotate(
    'text',
    x=50,
    y = .06,
    label='Nominal',
    color='red',
    family='Palatino Linotype'
  ) +
  annotate(
    'text',
    x=20,
    y=.055,
    label='Effective',
    color='blue',
    family='Palatino Linotype'
  )

p3 = ggplot() +
  geom_density(aes(x),color='red',fill='red',alpha=.5) +
  geom_density(aes(sjstats::weight(x,w_rfa)),
               color='blue',fill='blue',alpha=.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        text=element_text(family='Palatino Linotype')) +
  labs(
    x='Covariate',
    y='Density',
    title='RFA'
  ) +
  annotate(
    'text',
    x=50,
    y = .06,
    label='Nominal',
    color='red',
    family='Palatino Linotype'
  ) +
  annotate(
    'text',
    x=20,
    y=.055,
    label='Effective',
    color='blue',
    family='Palatino Linotype'
  )

grid.arrange(p1,p2,p3,
             layout_matrix=rbind(c(1,1,2,2),
                                 c(NA,3,3,NA)))



# monte carlo -------------------------------------------------------------

registerDoParallel(cores=4)
its = 500
set.seed(111)
foreach(i = 1:its,.combine='rbind') %dopar% {
  require(dplyr)
  require(randomForest)
  require(optmatch)
  require(lmtest)
  require(sandwich)
  n = 500 # n of observations
  ate = 5 # size of treatment effect
  
  # the confounder x
  x = rnorm(n = n, mean = 50, sd = 10)
  
  # the causal variable z
  prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
  z = rbinom(n = n, size = 1, prob = prZ)
  
  # the outcome variable y
  y = 1 + ate*z + 0.5*x + x^2 + rnorm(n = n, sd = 10)
  
  ndata = data.frame(y,z,x)
  
  # ols
  ols = lm(y~z*I(x - mean(x)),ndata)
  
  # matching
  mtc = lm(y~z+fullmatch(match_on(z~x),data=ndata),
           ndata)
  # rfa
  rdata = data.frame(
    y = ndata$y - suppressWarnings(predict(
      randomForest(y~x,data=ndata)
    )),
    z = ndata$z - suppressWarnings(predict(
      randomForest(z~x,data=ndata)
    ))
  )
  rfm = robustbase::lmrob(y~z,rdata)
  
  ros = robustbase::lmrob(rdata$y ~ z, ndata)
  
  # results
  return(rbind(
    coeftest(ols,vcovHC(ols)) %>%
      broom::tidy() %>%
      filter(term=='z') %>%
      mutate(model='Multiple Regression'),
    coeftest(mtc,vcovHC(mtc)) %>%
      broom::tidy() %>%
      filter(term=='z') %>%
      mutate(model='Matching'),
    coeftest(rfm) %>% #,vcovHC(rfm)) %>%
      broom::tidy() %>%
      filter(term=='z') %>%
      mutate(model='RFA'),
    coeftest(ros) %>%
      broom::tidy() %>%
      filter(term=='z') %>%
      mutate(model='Rosenbaum')
  ))
  #cat('progress:',round(100*i/its,2),'% \r')
} -> mc_rslts

mc_rslts %>%
  mutate(lo=estimate-1.96*std.error,
         hi=estimate+1.96*std.error) %>%
  group_by(model) %>%
  summarize(
    Bias = mean(estimate-ate),
    MSE = mean((estimate-ate)^2),
    Coverage = mean(lo<=ate & hi>=ate)
  ) %>%
  mutate_if(is.numeric,function(x)round(x,3)) %>%
  stargazer::stargazer(
    summary=F,
    rownames=F,
    header=F,
    title='Performance of Adjustment Strategies'
  )
mc_rslts %>%
  mutate(
    model=factor(
      model,
      levels=c('Multiple Regression','Matching','RFA','Rosenbaum')
    )
  ) %>%
  ggplot() +
  aes(estimate,model) +
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = ate) +
  theme_ridges() + 
  labs(
    x=expression(beta*' estimates'),
    y='',
    title='Distribution of ATEs'
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title.position = 'plot',
    text=element_text(family='Palatino Linotype')
  ) + 
  ggsave(
    'ate_distributions.png',
    height=3,
    width=6
  )

set.seed(111)
foreach(i = 1:its,.combine='rbind') %dopar% {
  require(dplyr)
  require(randomForest)
  require(optmatch)
  require(lmtest)
  require(sandwich)
  n = 500 # n of observations
  ate = 5 # size of treatment effect
  
  # the confounder x
  x = rnorm(n = n, mean = 50, sd = 10)
  
  # the causal variable z
  prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
  z = rbinom(n = n, size = 1, prob = prZ)
  
  # the outcome variable y
  y = 1 + ate*z + 0.5*x + x^2 + rnorm(n = n, sd = 10)
  
  # the data
  ndata = data.frame(y,z,x)
  
  # weights
  matches = fullmatch(z~x,data=ndata)
  w_ols = resid(lm(z~x + I(z*(x-mean(x))),ndata))^2
  w_mtc = resid(lm(z~matches,ndata))^2
  w_rfa = (
    z - predict(randomForest(z~x,data=ndata))
  )^2
  
  
  # results
  ovl_ols = overlap(x=x,w=w_ols)
  ovl_mtc = overlap(x=x,w=w_mtc)
  ovl_rfa = overlap(x=x,w=w_rfa)
  tibble(
    estimator = c('Multiple Regression','Matching','RFA'),
    '% overlap' = c(ovl_ols,ovl_mtc,ovl_rfa)*100
  )
  #cat('progress:',round(100*i/its,2),'% \r')
} -> mc_cvrg
mc_cvrg %>%
  mutate(
    estimator = factor(
      estimator,
      levels=c('Multiple Regression','Matching','RFA')
    )
  ) %>%
  ggplot() +
  aes(
    `% overlap`,
    estimator
  ) +
  geom_density_ridges(scale=1) +
  labs(
    y='',
    x='% Overlap',
    title='Overlap in Nominal & Effective Sample Distributions'
  ) +
  scale_x_continuous(breaks=seq(0,100,by=5)) +
  theme_ridges() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title.position = 'plot',
    text=element_text(family='Palatino Linotype')
  ) +
  ggsave(
    'overlap_distribution.png',
    height=3,
    width=6
  )
mc_cvrg %>%
  group_by(estimator) %>%
  summarize(mean=mean(`% overlap`))

set.seed(111)
foreach(i = 1:its,.combine='rbind') %dopar% {
  require(dplyr)
  require(randomForest)
  require(optmatch)
  require(lmtest)
  require(sandwich)
  n = 500 # n of observations
  ate = 5 # size of treatment effect
  
  # the confounder x
  x = rnorm(n = n, mean = 50, sd = 10)
  
  # the causal variable z
  prZ = 1/(1 + exp(2-(.05*(1.15*x - mean(x)))^2))
  z = rbinom(n = n, size = 1, prob = prZ)
  
  # the outcome variable y
  y = 1 + ate*z + 0.5*x + x^2 + rnorm(n = n, sd = 10)
  
  # the data
  ndata = data.frame(y,z,x)
  
  # weights
  matches = fullmatch(z~x,data=ndata)
  w_ols = resid(lm(z~x + I(z*(x-mean(x))),ndata))^2
  w_mtc = resid(lm(z~matches,ndata))^2
  w_rfa = (
    z - predict(randomForest(z~x,data=ndata))
  )^2
  
  
  # results
  ovl_ols = sample.mean(x=x,w=w_ols)
  ovl_mtc = sample.mean(x=x,w=w_mtc)
  ovl_rfa = sample.mean(x=x,w=w_rfa)
  tibble(
    method = c(
      rep('Regression',2),
      rep('Matching',2),
      rep('RFA',2)
    ),
    sample = c(
      rep(c('nominal','effective'),len=6)
    ),
    mean = c(
      mean_ols$sample_means[1],
      mean_ols$sample_means[2],
      mean_mtc$sample_means[1],
      mean_mtc$sample_means[2],
      mean_rfa$sample_means[1],
      mean_rfa$sample_means[2]
    )
  )
} -> mc_mean

mc_mean %>%
  mutate(
    its = rep(1:(n()/2),each=2)
  ) %>%
  group_by(its,method) %>%
  summarize(
    diff = diff(mean)
  ) %>%
  ggplot() +
  aes(
    diff,
    reorder(method,-diff)
  ) +
  geom_density_ridges(scale=1) +
  labs(
    y='',
    x='Difference in Effective and Nominal Sample Means',
    title='Effective vs. Nominal Sample Means'
  ) +
  theme_ridges() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title.position = 'plot',
    text=element_text(family='Palatino Linotype')
  ) +
  ggsave(
    'sample_means.png',
    height=3,
    width=6
  )



yres = y - 
  predict(
    randomForest(y~x,data=data)
  )
zres = z - predict(
    randomForest(z~x,data=data)
  )
summary(lm(yres~I(z-zhat)-1))

dd = sum(yres*zres)/sum(zres^2)
times(9999) %dopar% {
  obs = sample(1:length(yres))
  sum(yres*zres[obs])/sum(zres[obs]^2)
} -> permdd

ggplot() +
  geom_density(
    aes(permdd),
    col='black',
    fill='grey40'
  ) +
  geom_vline(
    xintercept=dd,col='red'
  )
test.ate = function(
  response,
  treatment,
  perms = 10000
) {
  
  # response
  y = response
  
  # treatment
  z = treatment
  
  # Estimate ATE
  d = sum(y*z)/sum(z^2)
  
  # Permute ATE
  times(perms) %dopar% {
    mix = sample(1:length(z))
    sum(y*z[mix])/sum(z[mix]^2)
  } -> permd
  
  # p.value
  p = 2*min(mean(d>=permd),mean(d<=permd))
  
  # Print Plot
  print(ggplot() +
    geom_density(
      aes(permd),
      col='black',
      fill='grey40'
    ) +
    geom_vline(
      xintercept=c(0,d),
      col=c('red','blue'),
      size=c(1,1)
    ) +
    labs(
      x = 'Possible Treatment Effects\nunder the Null',
      y = 'Probability',
      title = 'Null',
      subtitle = 'Observed'
    ) + theme_minimal() +
    theme(
      plot.title = element_text(
        size=12,color='red',
        hjust=(0-min(permd))/max(permd-min(permd))
      ),
      plot.subtitle = element_text(
        size=12,
        color='blue',
        hjust=(d-min(permd))/max(permd-min(permd))
      )
    ))
  
  # return results
  cat('     RESULTS FROM PERMUTATION TEST\n')
  cat('\n')
  cat('  Estimated Treatment Effect: ',round(d,3),'\n')
  cat('  Chance Probability:         ',round(p,3),'\n')
}
test.ate(
  response = yres,
  treatment = zres
)

cov(yres,zres)/cov(zres,zres)
sum((zres*yres)/(zres^2))
