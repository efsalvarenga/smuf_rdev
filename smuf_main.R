#===========================================
# Smart Metering Uncertainty Forecasting
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# smuf_main combined functions optim & fcst
#===========================================

#===========================================
# Libraries, Configs, Inputs
#===========================================
# library(forecast)
# library(verification)
# library(doParallel)
# library(rgenoud)
# library(rugarch)
# 
# cl  <- makeCluster(detectCores())
# registerDoParallel(cl)
# setwd("~/GitRepos/smuf_rdev")            # set working directory
# source("smuf_fxs.R")                     # Load functions script
# 
# wm01_00       <- readRDS("smuf_import-complete.rds")
# importpar     <- readRDS("smuf_import-parameter.rds")
# s01           <- importpar[1]
# s02           <- importpar[2]
# s03           <- importpar[3]
# sum_of_h      <- importpar[4]
# data_size     <- importpar[5]
# h             <- 0

#===========================================
# Integrated Parameters
#===========================================
cus_list      <- seq(1,30)
frontierstp   <- 5                       # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(0,1)*2069
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(2,3)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)


#===========================================
# BIG [h] LOOP Start
#===========================================
bighlpcrps = list()
bighlpsdev = list()

for (h in hrz_lim){
  ptm <- proc.time()
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  
  #===========================================
  # Individual customers forecast
  #===========================================
  cat("[Ind] ")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06       <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T,armalags,cross_overh)
  # wl06ag     <- fx_int_fcstgeneric_armagarch(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh)
  wv45       <- rowMeans(wl06[[1]])
  sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
  
  #===========================================
  # Random groups & evaluation
  #===========================================
  cat("[Rnd] ")
  wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
  wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
  wl06rnd    <- fx_int_fcst_kdcv(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T,armalags,cross_overh)
  wv45rnd    <- as.numeric(rowMeans(wl06rnd[[1]]) * rowSums(wm01_02l[[2]]))
  sd01rnd    <- as.numeric(fx_sd_mymat(wl06rnd[[3]]))
  cr01rnd    <- rowMeans(wl06rnd[[2]])
  
  #===========================================
  # Optimising groups & evaluation
  #===========================================
  wv46         <- seq(0,frontierstp)^2/frontierstp^2 * sum(wv45)
  cat("[OptSDEV] ")
  optgrp_sdev  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_sdev, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  res_sdev_kd  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcst_kdcv,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  res_sdev_ag  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  cat("[OptCVKD] ")
  optgrp_cvkd  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud","foreach"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_crps_kd, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  res_crps_kd  <- fx_applgrp(optgrp_cvkd,wv46,wm01_01,fx_int_fcst_kdcv,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  
  
  cat("[OptCVKD] ")
  optgrp_cvag  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud","foreach"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_crps_kd, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  res_crps_kd  <- fx_applgrp(optgrp_cvkd,wv46,wm01_01,fx_int_fcst_kdcv,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  
  
  
  
  
  
  
  
  
  
  bighlpcrps[[match(h,hrz_lim)]] = list(h,cbind(cr01rnd,wv45rnd),cbind(cr01optsdev,wv45optsdev),cbind(cr01optcrps,wv45optcrps))
  bighlpsdev[[match(h,hrz_lim)]] = list(h,cbind(sd01rnd,wv45rnd),cbind(sd01optsdev,wv45optsdev),cbind(sd01optcrps,wv45optcrps))
  print(proc.time() - ptm)
}
#===========================================
# BIG [h] LOOP Plots
#===========================================
for (i in 1:length(hrz_lim)){
  fx_plt_rnd_vs_opt(bighlpcrps[[i]],c(0.02,0.07),c(0,3),"CRPS")
}

#===========================================
# Outputs
#===========================================
# saveRDS(crpsh_CLU_dyn,  file="0200_analysis.rds") # crps winsize vs time_adead, mean&desv,                per grouping 
# saveRDS(wm01_4D[[1]],   file="0200_origdat1.rds") # data customer vs time_series,                         per grouping
# saveRDS(wm02_4D[[1]],   file="0200_seasona1.rds") # seas customer vs seasonal_ts,            per horizon, per grouping
# saveRDS(fcst_4D[[1]],   file="0200_forecas1.rds") # fcst winsize vs time_ahead per customer, per horizon, per grouping
# saveRDS(parbundl0200,   file="0200_parbundl.rds")
# saveRDS(kdscrps,        file="0200_function.rds")
