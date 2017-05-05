#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM0300 group optimisation
#-------------------------------------------
# Notes
# 11/04 Modular implemented
# 18/04 Parallel processing and data export
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
library(doParallel)
library(rgenoud)

ptm <- proc.time() # Start the clock!
cl  <- makeCluster(detectCores(),outfile="")
registerDoParallel(cl)

source('./KSSM0211_read-basefcst_v01.R')

#===========================================
# Optimisation functions
#===========================================
cus_grouping_A <- function(wv42){
  if (sum(wv42*wv45) > opt_min_cusd & sum(wv42*wv45) <= opt_max_cusd){
    fv01   <- (wv42 %*% (wm01_4D1[,(in_sample_ini):(event_horizon)])) / sum(wv42)
    fv02   <- decompose(msts(fv01[1,],seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
    fv03   <- fv02$random
    randev <- sd(fv03,na.rm=TRUE)
  } else {randev <- 10}#abs(sum(wv42*wv45))+1}
  return(randev)
}

cus_grouping_C <- function(wv42){
  if (sum(wv42*wv45) > opt_min_cusd & sum(wv42*wv45) <= opt_max_cusd){
    fl11   <- list((wv42 %*% (wm01_4D1[,(in_sample_ini):(event_horizon)])) / sum(wv42))
    fl12   <- kdscrps(fl11, c(0), win_size, ahead_t, (1 - crossvalsize / length(fv01) * s02), crossvalsize, sampling)
    fl13   <- ws_crossover_fx(fl12)
    randev <- mean(fl13[[1]][[1]][3,])
  } else {randev <- 10}#
  return(randev)
}

#===========================================
# Optimisation Parameters & Conf
#===========================================
opt_ahead_t   = 1
opt_win_sel   = 1
opt_hrz_sel   = 1
frontierstp   = detectCores() * 2

# Parameter Bundle
optgrppar     = c(opt_ahead_t,opt_win_sel,opt_hrz_sel,frontierstp)

# Optimisation Config
event_horizon = data_size*in_sample_fr+1 + hrz_lim[opt_hrz_sel]
in_sample_ini = event_horizon - event_horizon %/% s02 * s02 + 1
outsample_end = (data_size*(1-in_sample_fr)) %/% s02 * s02 + event_horizon - s02 - (hrz_lim[opt_hrz_sel]%/%s02*s02)
in_sample_siz = event_horizon - in_sample_ini + 1
outsample_siz = outsample_end - event_horizon
inosample_siz = in_sample_siz + outsample_siz
event_hrz_mod = event_horizon - in_sample_ini + 1

wv42          = cus_list * 0
wv43          = cus_list * 0
for (j in cus_list){
  wv43[j] = fcst_4D1[[opt_hrz_sel]][[j]][opt_win_sel,opt_ahead_t]
}
wv44          = wm02_4D1[[opt_hrz_sel]][,opt_ahead_t]
wv45          = wv43+wv44                          # demand forecast for opt selection
wv46          = seq(0,sum(wv45),sum(wv45)/frontierstp)

#===========================================
# Optimisation Run
#===========================================
optgrp_plA <- foreach (i = 1:frontierstp,
                       .packages=c("forecast","rgenoud"),
                       .combine=c("rbind")) %dopar%{
  opt_min_cusd  = wv46[i]
  opt_max_cusd  = wv46[i+1]
  optgrp   <- genoud(cus_grouping_A, nvars=length(cus_list), max.generations=300, wait.generations=20,
                     starting.values=c(rep(1,length(cus_list))),
                     Domains = cbind(c(rep(0,length(cus_list))),c(rep(1,length(cus_list)))),
                     data.type.int=TRUE,  int.seed=1,
                     print.level=0)
  cat("\nDone stepwise frontier",i,"of",frontierstp)
  optgrp$par
}
optgrp_plA = matrix(optgrp_plA,frontierstp,length(cus_list))

optgrp_plC <- foreach (i = 1:frontierstp,
                       .packages=c("forecast","verification","rgenoud"),
                       .combine=c("rbind")) %:%{
                         opt_min_cusd  = wv46[i]
                         opt_max_cusd  = wv46[i+1]
                         optgrp   <- genoud(cus_grouping_C, nvars=length(cus_list), max.generations=300, wait.generations=20,
                                            starting.values=c(rep(1,length(cus_list))),
                                            Domains = cbind(c(rep(0,length(cus_list))),c(rep(1,length(cus_list)))),
                                            data.type.int=TRUE,  int.seed=1,
                                            print.level=1)
                         cat("\nDone stepwise frontier",i,"of",frontierstp)
                         optgrp$par
                       }
optgrp_plC = matrix(optgrp_plC,frontierstp,length(cus_list))

# optgrp_plt = matrix(nrow=frontierstp,ncol=3)
# for (i in 1:frontierstp){
#   opt_min_cusd  = wv46[i]
#   opt_max_cusd  = wv46[i+1]
#   optgrp_plt[i,1] = sum(optgrp_plA[i,]*wv45)
#   optgrp_plt[i,2] = cus_grouping_A(optgrp_plA[i,])
#   optgrp_plt[i,3] = sum(optgrp_plA[i,])
# }

print(proc.time() - ptm)        # Stop the clock

#===========================================
# Outputs
#===========================================
saveRDS(list(optgrp_plA,optgrp_plt), file="0300_optgrp.rds")
saveRDS(optgrppar,                   file="0300_optpar.rds")