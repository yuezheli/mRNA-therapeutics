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
mod <- mread("varga2005_m1") %>% init(Complex_extracellular = 5e4)

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


set1 <- c('k_escape')
set2 <- c('k_unpack', 'k_NPC','k_in' ,'k_dissociation')

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
    summarise(pexp = auc_partial(time, Plasmid_nuclear)) %>% 
    pull(pexp)
}

##------------------------- Global sensitivity analysis -------------------------##
# warning: the following line takes several minutes to run
pglobal = sobol2007(batch_mRNA, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot)
plot(pglobal)
lines(c(0,7), c(0.05, 0.05), type="l", lty = 2)

