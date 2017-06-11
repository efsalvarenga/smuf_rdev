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
# Initialising
#===========================================
setwd("~/GitRepos/smuf_rdev")
source("smuf_fxs.R")
savfile = "smuf_main_XXXX.rds"

# From smuf_import
wm01_00       <- readRDS("smuf_import-complete.rds")
importpar     <- readRDS("smuf_import-parameter.rds")
s01           <- importpar[1]
s02           <- importpar[2]
s03           <- importpar[3]
sum_of_h      <- importpar[4]
data_size     <- importpar[5]

#===========================================
# Integrated Parameters
#===========================================
cus_list      <- seq(1,100)
frontierstp   <- 20                      # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- 0 #seq(0,1)*2069
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 16                      # Steps used for multiple crossvalidation (Only KDE)
crossvalfocus <- max(ahead_t)            # What period is focused when running crossvalidation
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(8,8)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)


#===========================================
# BIG [h] LOOP Start
#===========================================
bighlpopgr <- list()
bighlpcrps <- list()

for (h in hrz_lim){
  ptm    <- proc.time()
  runkey <- Sys.time()
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  
  #===========================================
  # Individual customers forecast
  #===========================================
  cat("[Ind] ")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06       <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T,armalags,cross_overh)
  wv45       <- rowMeans(wl06[[1]])
  sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
  wv46       <- seq(0,frontierstp)^2/frontierstp^2 * sum(wv45)
  
  #===========================================
  # Random groups & evaluation
  #===========================================
  cat("[Rnd] ")
  wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
  wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
  wl06rnd    <- fx_int_fcst_kdcv(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T,armalags,cross_overh)
  wv45rnd    <- as.numeric(rowMeans(wl06rnd[[1]]) * rowSums(wm01_02l[[2]]))
  # sd01rnd    <- as.numeric(fx_sd_mymat(wl06rnd[[3]]))
  cr01rnd    <- rowMeans(wl06rnd[[2]])
  
  #===========================================
  # Optimising groups & evaluation
  #===========================================
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
  bighlpopgr   <- fx_sav_optgrps(c("sdev",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_sdev)
  res_sdev_kd  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcst_kdcv,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  res_sdev_ag  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  
  cat("[OptCVKD] ")
  optgrp_cvkd  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud","foreach"),
                           .export=c('fcst_mccv'),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_crps, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1, print.level=1,
                                                use_arma=F)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  bighlpopgr   <- fx_sav_optgrps(c("cvkd",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_cvkd)
  res_crps_kd  <- fx_applgrp(optgrp_cvkd,wv46,wm01_01,fx_int_fcst_kdcv,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  res_crps_ag  <- fx_applgrp(optgrp_cvkd,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  
  # cat("[OptCVAG] ")
  # optgrp_cvag  <- foreach (i = 1:frontierstp,
  #                          .packages=c("forecast","rgenoud","foreach"),
  #                          .export=c('fcst_mccv'),
  #                          .combine=c("rbind")) %dopar% {
  #                            opt_min_cusd  = wv46[i]
  #                            opt_max_cusd  = wv46[i+1]
  #                            optgrp   <- genoud(fx_optgrp_crps, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
  #                                               starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
  #                                               data.type.int=TRUE,  int.seed=1, print.level=1,
  #                                               use_arma=T)
  #                            if(optgrp$value == 10) {
  #                              grouped = c(rep(0,nrow(wm01_01)))
  #                            } else {
  #                              grouped = optgrp$par
  #                            }
  #                            grouped
  #                          }
  # bighlpopgr   <- fx_sav_optgrps(c("cvag",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_cvag)
  # res_crps_ag  <- fx_applgrp(optgrp_cvag,wv46,wm01_01,fx_int_fcst_kdcv,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)

  bighlpcrps   <- fx_sav_optress(c("sdev+crps_kd+ag",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),
                                 list(c(h,frontierstp,length(cus_list)),cbind(cr01rnd,wv45rnd),res_sdev_kd,res_sdev_ag,res_crps_kd,res_crps_ag))
  fx_plt_rnd_vs_opt(bighlpcrps[[length(bighlpcrps)]][[2]],c(0.01,0.04),c(0,8),"CRPS")
  print(proc.time() - ptm)
}

#===========================================
# Outputs
#===========================================
saveRDS(list(bighlpopgr,bighlpcrps),  file=savfile)