#===========================================
# Smart Metering Uncertainty Forecasting
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# smuf_lags checking lags for ARMA-GARCH
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
library(forecast)
library(doParallel)
library(rugarch)

ptm <- proc.time() # Start the clock!
cl  <- makeCluster(detectCores())
registerDoParallel(cl)
setwd("~/GitRepos/smuf_rdev")

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
cus_list      <- seq(1,10)
# frontierstp   <- 5                       # Number of demand bins (Stepwise frontier for portfolio optimisation)
# win_size      <- c(4,24)                 # Small and large win_size (select only 2)
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(0,1)*2537
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
# crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
# crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
seas_bloc_ws  <- 6                       # Number of weeks used for calculating seasonality pattern (6 seems best)
# sampling      <- 1024                    # For monte-carlo CRPS calculation
maxlag        <- 5                       # Max lags analysed for ARIMA fit (ARMA-GARCH model)

#===========================================
# Functions Declarations: Modules
#===========================================
fx_lags_armagarch <- function (wm04,maxlag){#,ahead_t,out_evhor,sampling){
  lags_armagarch <- foreach (j = 1:nrow(wm04), .packages=c("rugarch")) %dopar% {
    runvec       <- wm04[j,1:out_evhor[7]]
    # Defining ARMA lags
    final.bic <- matrix(nrow=0,ncol=4)
    for (p in 0:maxlag) for (q in 0:maxlag) {
      if ( p == 0 && q == 0) {
        next
      }
      arimaFit = tryCatch(arima(runvec, order=c(p, 0, q)),
                          error=function( err ) FALSE,
                          warning=function( err ) FALSE )
      if( !is.logical( arimaFit ) ) {
        final.bic <- rbind(final.bic,c(p,q,BIC(arimaFit),AIC(arimaFit)))
      } else {
        next
      }
    }
    final.ord <- final.bic[sort.list(final.bic[,3]), ]
  }
  return(fcst_armagarch)
}

#===========================================
# BIG [h] LOOP Start
#===========================================
wm01_01    <- wm01_00[min(cus_list):length(cus_list),]
out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,seas_bloc_ws,0)
wm01       <- wm01_01[,out_evhor[2]:out_evhor[3]]                       # work matrix
wm02       <- fx_seas(wm01,s01,s02,sum_of_h,out_evhor)                  # in-sample seasonality pattern
wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]                      # out-sample original load data
wm04       <- fx_unseas(wm01,wm02,s02,out_evhor)                        # in-out sample unseasonalised



bighlpcrps = list()
bighlpsdev = list()

for (h in hrz_lim){
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  #===========================================
  # Individual customers forecast
  #===========================================
  cat("[Ind] ")
  wm01_01    <- wm01_00[min(cus_list):length(cus_list),]
  wl06       <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,seas_bloc_ws,crossvalsize,T)
  wv45       <- rowMeans(wl06[[1]])
  sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
  
  #===========================================
  # Random groups forecast
  #===========================================
  cat("[Rnd] ")
  wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
  wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
  wl06rnd    <- fx_int_fcst_kdcv(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_size,seas_bloc_ws,crossvalsize,T)
  wv45rnd    <- as.numeric(rowMeans(wl06rnd[[1]]) * rowSums(wm01_02l[[2]]))
  sd01rnd    <- as.numeric(fx_sd_mymat(wl06rnd[[3]]))
  cr01rnd    <- rowMeans(wl06rnd[[2]])
  
  #===========================================
  # Optimised sdev groups forecast
  #===========================================
  cat("[OptSDEV] ")
  wv46         <- seq(0,frontierstp)^2/frontierstp^2 * sum(wv45)
  optgrp_sdev  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_sdev, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             optgrp$par
                           }
  opt_min_cusd<- 0
  opt_max_cusd<- max(wv46)
  wm01_03l    <- list(optgrp_sdev %*% wm01_01, optgrp_sdev)
  wm01_03     <- wm01_03l[[1]] / rowSums(wm01_03l[[2]])
  wl06optsdev <- fx_int_fcst_kdcv(wm01_03,h,in_sample_fr,s01,s02,sum_of_h,win_size,seas_bloc_ws,crossvalsize,T)
  wv45optsdev <- as.numeric(rowMeans(wl06optsdev[[1]]) * rowSums(wm01_03l[[2]]))
  sd01optsdev <- as.numeric(fx_sd_mymat(wl06optsdev[[1]]))
  cr01optsdev <- rowMeans(wl06optsdev[[2]])
  cr02optsdev <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud","foreach"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_crps(optgrp_sdev[i,])
                          }
  cr02optsdev <- as.numeric(cr02optsdev)
  sd02optsdev <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_sdev(optgrp_sdev[i,])
                          }
  sd02optsdev <- as.numeric(sd02optsdev)
  
  # for sdev calc: sd01 is the outsample result, sd02 is the crossval result
  # for crps calc: cr01 is the outsample result, cr02 is the crossval result
  
  #===========================================
  # Optimised crps groups forecast
  #===========================================
  cat("[OptCrossval]\n")
  optgrp_crps  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud","foreach"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_crps, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             optgrp$par
                           }
  opt_min_cusd<- 0
  opt_max_cusd<- max(wv46)
  wm01_04l    <- list(optgrp_crps %*% wm01_01, optgrp_crps)
  wm01_04     <- wm01_04l[[1]] / rowSums(wm01_04l[[2]])
  wl06optcrps <- fx_int_fcst_kdcv(wm01_04,h,in_sample_fr,s01,s02,sum_of_h,win_size,seas_bloc_ws,crossvalsize,T)
  wv45optcrps <- as.numeric(rowMeans(wl06optcrps[[1]]) * rowSums(wm01_04l[[2]]))
  sd01optcrps <- as.numeric(fx_sd_mymat(wl06optcrps[[1]]))
  cr01optcrps <- rowMeans(wl06optcrps[[2]])
  cr02optcrps <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud","foreach"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_crps(optgrp_crps[i,])
                          }
  cr02optcrps <- as.numeric(cr02optcrps)
  sd02optcrps <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_sdev(optgrp_crps[i,])
                          }
  sd02optcrps <- as.numeric(sd02optcrps)
  
  # for sdev calc: sd01 is the outsample result, sd02 is the crossval result
  # for crps calc: cr01 is the outsample result, cr02 is the crossval result
  
  bighlpcrps[[match(h,hrz_lim)]] = list(h,cbind(cr01rnd,wv45rnd),cbind(cr01optsdev,wv45optsdev),cbind(cr02optsdev,wv45optsdev),cbind(cr01optcrps,wv45optcrps),cbind(cr02optcrps,wv45optcrps))
  bighlpsdev[[match(h,hrz_lim)]] = list(h,cbind(sd01rnd,wv45rnd),cbind(sd01optsdev,wv45optsdev),cbind(sd02optsdev,wv45optsdev),cbind(sd01optcrps,wv45optcrps),cbind(sd02optcrps,wv45optcrps))
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
