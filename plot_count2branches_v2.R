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

##reorder factor for plotting
data$taxon3 <- factor(data$taxon3, levels=c('sm', 'lg'))

##get models
##set the theme so it's not hideous
theme_set(theme_tidybayes() + panel_border())

##using tidybayes
##see https://mjskay.github.io/tidybayes/articles/tidybayes.html
data %>%
  group_by(taxon) %>%
  data_grid(blsr1 = seq_range(blsr1, n = 228), taxon2 = taxon2, taxon3 = taxon3, lab = lab) %>%
  add_epred_draws(fit2) %>%
  ggplot(aes(x = blsr1, y = count, color = taxon)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95)) +
  geom_point(data = data, aes(color = taxon2, size = taxon3)) +
  geom_text(data = data, na.rm=T,  nudge_x=0.01, nudge_y=-2, aes(label = lab)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_discrete(type = c("#0000ffff", "#d40000ff", "grey30", "#00000087"), breaks = c("ancestral", "bat", "not")) + 
  theme(legend.position = c(0.8, 0.2), legend.title=element_blank()) + 
  guides(size="none", fill="none", text="none") + 
  labs(x='Branch Lengths calculated from 4d sites', y='No. genes under selection in GO:0002376')
  
##save
ggsave("count2branches_phylofit.pdf", w=7, h=5)

data %>%
  group_by(taxon) %>%
  data_grid(blsr2 = seq_range(blsr2, n = 228), taxon2 = taxon2, taxon3 = taxon3, lab = lab) %>%
  add_epred_draws(fit2i) %>%
  ggplot(aes(x = blsr2, y = count, color = taxon)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95)) +
  geom_point(data = data, aes(color = taxon2, size = taxon3)) +
  geom_text(data = data, na.rm=T,  nudge_x=0.01, nudge_y=-2, aes(label = lab)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_discrete(type = c("#0000ffff", "#d40000ff", "grey30", "#00000087"), breaks = c("ancestral", "bat", "not")) + 
  theme(legend.position = c(0.8, 0.2), legend.title=element_blank()) + 
  guides(size="none", fill="none", text="none") + 
  labs(x='Branch Lengths (coalescent units)', y='No. genes under selection in GO:0002376')
  
##save
ggsave("count2branches_iqTree.pdf", w=7, h=5)

data %>%
  group_by(taxon) %>%
  data_grid(blma = seq_range(blma, n = 228), taxon2 = taxon2, taxon3 = taxon3, lab = lab) %>%
  add_epred_draws(fit2m) %>%
  ggplot(aes(x = blma, y = count, color = taxon)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95)) +
  geom_point(data = data, aes(color = taxon2, size = taxon3)) +
  geom_text(data = data, na.rm=T,  nudge_x=3, nudge_y=-2, aes(label = lab)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_discrete(type = c("#0000ffff", "#d40000ff", "grey30", "#00000087"), breaks = c("ancestral", "bat", "not")) + 
  theme(legend.position = c(0.8, 0.2), legend.title=element_blank()) + 
  guides(size="none", fill="none", text="none") + 
  labs(x='Branch Lengths (mya)', y='No. genes under selection in GO:0002376')
  
##save
ggsave("count2branches_Ma_simple.pdf", w=7, h=5)

data %>%
  group_by(taxon) %>%
  data_grid(blma = seq_range(blma, n = 228), taxon2 = taxon2, taxon3 = taxon3, lab = lab) %>%
  add_epred_draws(fit4m) %>%
  ggplot(aes(x = blma, y = count, color = taxon)) +
  stat_lineribbon(aes(y = .epred), .width = c(0.95)) +
  geom_point(data = data, aes(color = taxon2, size = taxon3)) +
  geom_text(data = data, na.rm=T,  nudge_x=3, nudge_y=-2, aes(label = lab)) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_discrete(type = c("#0000ffff", "#d40000ff", "grey30", "#00000087"), breaks = c("ancestral", "bat", "not")) + 
  theme(legend.position = c(0.8, 0.2), legend.title=element_blank()) + 
  guides(size="none", fill="none", text="none") + 
  labs(x='Branch Lengths (mya)', y='No. genes under selection in GO:0002376')
  
##save
ggsave("count2branches_Ma_complex.pdf", w=7, h=5)

