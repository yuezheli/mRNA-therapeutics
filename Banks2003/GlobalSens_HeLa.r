# Global sensitivity analysis, HeLa cells
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(sensitivity)
library(PKPDmisc)
library(mrgsim.parallel)

mod <- mread("banks2003") %>% init(M = 2.41e11)

# set up sampling for global sensitivity analysis
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

set2 <- c('k1','k2','k3','k4')
N <- sampleperparam * length(set2)
sets <- list(
  list(which = set2, factor = c(0.1,10), N = N, l = param(mod))
)
samps <- map(sets, do.call, what = gen_samples)
samp <- map(c(1,2), ~ map_dfc(samps, .x))

batch_nPlasmid <- function(x) {
  mod %>% 
    idata_set(x) %>%
    mrgsim(obsonly = TRUE) %>% 
    group_by(ID) %>% 
    summarise(pexp = auc_partial(time, N)) %>% 
    pull(pexp)
}

##------------------------- Global sensitivity analysis -------------------------##
# warning: the following line takes several minutes to run
pglobal = sobol2007(batch_nPlasmid, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot)
plot(pglobal)
lines(c(0,7), c(0.05, 0.05), type="l", lty = 2)
# save a minimal set of the Sobol analysis data
# that is sufficient for visualization
sobolx <- list(S = pglobal$S , T = pglobal$T)
saveRDS(sobolx, file = "data/Sobol_nucleusPlasmid.rds")


##------------------------- Visualization -------------------------##

x <- readRDS("data/Sobol_nucleusPlasmid.rds")

set2 <- c('k1', 'k2', 'k3','k4')

globSens_nP <- tibble(Parameter = set2,
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

fig_sens_glob_nP <- ggplot(data=globSens_nP, aes(x=Parameter, y=Index, group=Effect, col=Effect)) +
  geom_point(position = position_dodge(width=0.3)) +
  geom_errorbar(aes(ymax=hi, ymin=lo), position=position_dodge(width=0.3), width=0) +
  theme_bw() +
  geom_hline(aes(yintercept=0.05), lty=2) +
  theme(legend.title = element_text(size = 15), legend.text=element_text(size=12)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12)) + 
  ggtitle('Global sensitivity analysis, use nucleus plasmid as a read out')

print(fig_sens_glob_nP)
ggsave(filename = 'img/GlobalSensHeLa.png', plot = fig_sens_glob_nP, width = 10, height = 3)
