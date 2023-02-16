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

##clear out cache
rm(list=ls())

## data
load("count2branches.RData")

##print out the comparison
sink("model_comparison_blsr_phylofit.txt")
print(loo_compare(fit0,fit.p,fit1, fit2, fit3, fit4, criterion = "waic"))
sink()

##print out best-fit model
sink("best_model_blsr_phylofit.txt")
print(summary(fit2))
print(bayes_R2(fit2))
sink()

##print out the comparison
sink("model_comparison_blsr_iqTree.txt")
print(loo_compare(fit0i, fit.pi, fit1i, fit2i, fit3i, fit4i, criterion = "waic"))
sink()

##print out best-fit model
sink("best_model_blsr_iqTree.txt")
print(summary(fit2i))
print(bayes_R2(fit2i))
sink()

##print out the comparison
sink("model_comparison_blsr_Ma.txt")
print(loo_compare(fit0m, fit.pm, fit1m, fit2m, fit3m, fit4m, criterion = "waic"))
sink()

##print out best-fit model
sink("best_model_blsr_Ma.txt")
print(summary(fit2m))
print(bayes_R2(fit2m))
sink()

##print out comparison across preditors
sink("model_comparison_predictors.txt")
print(loo_compare(fit2, fit2i, fit2m, fit4m, criterion = "waic"))
sink()

##all sources of variance accounted for
##https://discourse.mc-stan.org/t/difference-between-method-fitted-vs-predict-in-marginal-effects-function/6901/3
##assemble predictions for Ariadna

pred <- as.data.frame(cbind(data$label, data$count, predict(fit2, robust=T)[,1], predict(fit2i, robust=T)[,1], predict(fit2m, robust=T)[,1], predict(fit1, robust=T)[,1], predict(fit1i, robust=T)[,1], predict(fit1m, robust=T)[,1]))

##add names
colnames(pred) <- c("label", "count", 	"phylofit 2 intercept", "iqTree 2 intercept", "Ma 2 intercept", "phylofit 1 intercept", "iqTree 1 intercept", "Ma 1 intercept")

##export
write.csv(pred, "predicted.csv", row.names = F)
