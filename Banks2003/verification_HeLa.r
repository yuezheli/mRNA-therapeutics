rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at current folder

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

# read in observed data
obs_fig3 <- read.csv('data/Banks2003_fig3.csv', header = TRUE)

# simulation
mod <- mread("banks2003") %>% init(M = 2.41e11)

sim0 <- mod %>% mrgsim(delta = 0.1, end = 10) %>% as.tibble()

# plot result
mp <- ggplot(data = sim0, aes(x = time)) + 
  geom_line(aes(y = M, color = 'simul')) + 
  geom_point(data = obs_fig3 %>% filter(type == 'medium'), aes(y = value, color = 'Banks 2003')) + 
  xlim(0, 10) + ylim(2.406e11, 2.411e11) + scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  labs(x = 'time (h)', y = 'medium plasmid count', color = '') + theme_bw() +  theme(legend.position = "bottom")

cp <- ggplot(data = sim0, aes(x = time)) + 
  geom_line(aes(y = C, color = 'simul')) + 
  geom_point(data = obs_fig3 %>% filter(type == 'cytosol'), aes(y = value, color = 'Banks 2003')) + 
  xlim(0, 10) + ylim(0, 2e8) + scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  labs(x = 'time (h)', y = 'cytosol plasmid count', color = '') + theme_bw() +  theme(legend.position = "bottom")

np <- ggplot(data = sim0, aes(x = time)) + 
  geom_line(aes(y = N, color = 'simul')) + 
  geom_point(data = obs_fig3 %>% filter(type == 'nucleus'), aes(y = value, color = 'Banks 2003')) + 
  xlim(0, 10) + ylim(0, 2.5e8) + scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  labs(x = 'time (h)', y = 'cell nucleus plasmid count', color = '') + theme_bw() +  theme(legend.position = "bottom")

# save output
png('img/verification_HeLa.png', width = 900, height = 300)
grid.arrange(mp, cp, np, nrow = 1)
dev.off()