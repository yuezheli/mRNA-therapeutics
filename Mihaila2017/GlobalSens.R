rm(list = ls())
gc()
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script
# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(mrggsave)
library(sensitivity)
library(PKPDmisc)
library(mrgsim.parallel)
set.seed(88771) 

##------------------------- Prepare model -------------------------##
# add volume
Vextra = 3e-4 # extracellular compartment volume; unit L-1
Vintra = 1.4e-12 # intracellular compartment volume; unit L-1

# constant conversion
Avogadro = 6.02e23

# convert initial condition
InitConv <- function(MolCount = 10, # unit: molecule count; output: nM
                     Volume = Vintra, const = Avogadro){ return(MolCount/ Avogadro/ Vintra * 1e9) }

init3 <- list(
  E = 2.8e7,  # nM
  R = InitConv(1e4), # nM 
  M = InitConv(100) # nM
)

mod <- mread("mihaila2017_v3")  %>% init(init3) 

##------------------------- Prepare global sensitivity analysis -------------------------##

# create sampling method
gen_samples <- function(n, l, which = names(l), 
                        factor = c(0.01,100), N = NULL) { # NEW ARGUMENT: N, the absolute sampling number per parameter
  
  vars <- tidyselect::vars_select(names(l), !!(enquo(which)))
  
  l <- as.list(l)[vars]
  
  l <- lapply(l, function(x) x*factor )
  
  if(is.numeric(N)) { # NEW
    n <- N*2  
  } else {
    n <- length(l)*n*2
  }
  
  df <- as.data.frame(l)
  
  len <- length(df)
  
  X <- matrix(ncol=len, nrow=n)
  
  colnames(X) <- names(df)
  
  Y <- X
  
  for(i in seq(len)){
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    X[,i] <- r
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    Y[,i] <- r
  }
  
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}
# set up simulation parameters for sobol2007 function
simulationboot = 500 # number of simulation time in sobol2007
sampleperparam = 150 # sampling size per parameter, for sobol2007


set1 <- c('k1')
set2 <- c('k2', 'k3','k4' ,'k5', 'k6', 'k7')

N <- sampleperparam * length(c(set1,set2))
# set different variation in each group
sets <- list(
  list(which = set1, factor = c(0.5, 1.5), N = N, l = param(mod)), 
  list(which = set2, factor = c(0.01,100), N = N, l = param(mod))
)

samps <- map(sets, do.call, what = gen_samples)
samp <- map(c(1,2), ~ map_dfc(samps, .x))

batch_mRNA <- function(x) {
  mod %>% 
    idata_set(x) %>%
    mrgsim(obsonly = TRUE) %>% 
    group_by(ID) %>% 
    summarise(pexp = auc_partial(time, RNAcount)) %>% 
    pull(pexp)
}

##------------------------- Global sensitivity analysis -------------------------##
# warning: the following line takes several minutes to run
pglobal = sobol2007(batch_mRNA, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot)
plot(pglobal)
lines(c(0,7), c(0.05, 0.05), type="l", lty = 2)
# save a minimal set of the Sobol analysis data
# that is sufficient for visualization
sobolx <- list(S = pglobal$S , T = pglobal$T)
saveRDS(sobolx, file = "data/Sobol_mRNA.rds")

##------------------------- Visualization -------------------------##

x <- readRDS("data/Sobol_mRNA.rds")

set1 <- c('k1')
set2 <- c('k2', 'k3','k4' ,'k5', 'k6', 'k7')

globSens_mRNA <- tibble(Parameter = c(set1,set2),
                        main = x$S$original,
                        main_lo = x$S$'min. c.i.',
                        main_hi = x$S$'max. c.i.',
                        total = x$T$original,
                        total_lo = x$T$'min. c.i.',
                        total_hi = x$T$'max. c.i.') %>%
  gather(effect, Index, -Parameter, -main_lo, -main_hi, -total_lo, -total_hi) %>%
  mutate(lo = ifelse(effect == "main", main_lo, total_lo),
         hi = ifelse(effect == "main", main_hi, total_hi),
         Effect = factor(effect))

fig_sens_glob_mRNA <- ggplot(data=globSens_mRNA, aes(x=Parameter, y=Index, group=Effect, col=Effect)) +
  geom_point(position = position_dodge(width=0.3)) +
  geom_errorbar(aes(ymax=hi, ymin=lo), position=position_dodge(width=0.3), width=0) +
  theme_bw() +
  geom_hline(aes(yintercept=0.05), lty=2) +
  theme(legend.title = element_text(size = 15), legend.text=element_text(size=12)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12)) + 
  ggtitle('Global sensitivity analysis, use mRNA as a read out')

print(fig_sens_glob_mRNA)

