rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(tidyverse)
library(mrgsolve)

# list time points to save 

savetime <- seq(from = 0, to = 72, by = 1) * 3600

mod <- mread("../model2")

simul <- mod %>% mrgsim() %>% as.tibble() %>% 
         filter(time %in% savetime) %>% 
         mutate(hour = time/3600) %>%
         select(-c(ID, PlasmaDrug, mRNA, time))

# save the result
write.csv(simul, file = "../data/ForJuliaValidation.csv", row.names = FALSE)