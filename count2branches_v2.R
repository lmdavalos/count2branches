##modeling relationship PSGs to branch lengths

##run in stan
library(brms)

##use for plotting
library(magrittr)
library(dplyr)
library(forcats)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(emmeans)
library(broom)
library(rstan)
library(rstanarm)
library(bayesplot)
library(MCMCglmm)
library(RColorBrewer)

##get data lightly edited
data<-read.csv("summary.fittedValues_ed.csv")

##order factor
##helps make relative to non bats
data$taxon <- factor(data$taxon, levels=c('not', 'bat'))

##this series of models uses branch lengths from phylofit
##one intercept and slope
fit1 <- brm(bf(count ~ sqrt(blsr1)), data = data, family= negbinomial(link='log'), chains = 4, cores = 4)

##two intercepts, one slope--ancova
fit2 <- brm(bf(count ~ sqrt(blsr1) + taxon), data = data, family= negbinomial(link='log'),
 chains = 4, cores = 4)
 
##two intercepts, two slopes
fit3 <- brm(bf(count ~ sqrt(blsr1) : taxon), data = data, family= negbinomial(link='log'),
 chains = 4, cores = 4)

##one intercept, two slopes
fit4 <- brm(bf(count ~ sqrt(blsr1) * taxon), data = data, family= negbinomial(link='log'),  chains = 4, cores = 4)

##of these, fit2 is best so run those with linear and Poisson alternatives
##Poisson wih sqrt
fit.p <- brm(bf(count ~sqrt(blsr1) + taxon), data = data, family= poisson(link='log'),  chains = 4, cores = 4)

## neg binomial no sqrt transformation
fit0 <- brm(bf(count ~ blsr1 + taxon), data = data, family= negbinomial(link='log'),  chains = 4, cores = 4)

##add waic to all models
fit1 <- add_criterion(fit1, "waic")
fit2 <- add_criterion(fit2, "waic")
fit3 <- add_criterion(fit3, "waic")
fit4 <- add_criterion(fit4, "waic")
fit.p <- add_criterion(fit.p, "waic")
fit0 <- add_criterion(fit0, "waic")

##save
save.image("count2branches.RData")

##this series of models uses branch lengths from iqtree
##one intercept and slope
fit1i <- brm(bf(count ~ sqrt(blsr2)), data = data, family= negbinomial(link='log'), chains = 4, cores = 4)

##two intercepts, one slope--ancova
fit2i <- brm(bf(count ~ sqrt(blsr2) + taxon), data = data, family= negbinomial(link='log'),
 chains = 4, cores = 4)
 
##two intercepts, two slopes
fit3i <- brm(bf(count ~ sqrt(blsr2) : taxon), data = data, family= negbinomial(link='log'),
 chains = 4, cores = 4)

##one intercept, two slopes
fit4i <- brm(bf(count ~ sqrt(blsr2) * taxon), data = data, family= negbinomial(link='log'),  chains = 4, cores = 4)

##of these, fit2 is best so run those with linear and Poisson alternatives
##Poisson wih sqrt
fit.pi <- brm(bf(count ~sqrt(blsr2) + taxon), data = data, family= poisson(link='log'),  chains = 4, cores = 4)

## neg binomial no sqrt transformation
fit0i <- brm(bf(count ~ blsr2 + taxon), data = data, family= negbinomial(link='log'),  chains = 4, cores = 4)

##add waic to all models
fit1i <- add_criterion(fit1i, "waic")
fit2i <- add_criterion(fit2i, "waic")
fit3i <- add_criterion(fit3i, "waic")
fit4i <- add_criterion(fit4i, "waic")
fit.pi <- add_criterion(fit.pi, "waic")
fit0i <- add_criterion(fit0i, "waic")

##save
save.image("count2branches.RData")

##this series of models uses branch lengths from iqtree
##one intercept and slope
fit1m <- brm(bf(count ~ sqrt(blma)), data = data, family= negbinomial(link='log'), chains = 4, cores = 4)

##two intercepts, one slope--ancova
fit2m <- brm(bf(count ~ sqrt(blma) + taxon), data = data, family= negbinomial(link='log'),
 chains = 4, cores = 4)
 
##two intercepts, two slopes
fit3m <- brm(bf(count ~ sqrt(blma) : taxon), data = data, family= negbinomial(link='log'),
 chains = 4, cores = 4)

##one intercept, two slopes
fit4m <- brm(bf(count ~ sqrt(blma) * taxon), data = data, family= negbinomial(link='log'),  chains = 4, cores = 4)

##of these, fit2 is best so run those with linear and Poisson alternatives
##Poisson wih sqrt
fit.pm <- brm(bf(count ~sqrt(blma) + taxon), data = data, family= poisson(link='log'),  chains = 4, cores = 4)

## neg binomial no sqrt transformation
fit0m <- brm(bf(count ~ blma + taxon), data = data, family= negbinomial(link='log'),  chains = 4, cores = 4)

##add waic to all models
fit1m <- add_criterion(fit1m, "waic")
fit2m <- add_criterion(fit2m, "waic")
fit3m <- add_criterion(fit3m, "waic")
fit4m <- add_criterion(fit4m, "waic")
fit.pm <- add_criterion(fit.pm, "waic")
fit0m <- add_criterion(fit0m, "waic")

##save
save.image("count2branches.RData")
