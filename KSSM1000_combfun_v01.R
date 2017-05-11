#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM1000 combined functions optim & fcst
#-------------------------------------------
# Notes
# 11/04 Modular implemented
# 15/04 Saves forecasts and data for later use
# 26/04 Cross validation and seasonal bloc
# 10/05 Finished optimisation grouping
# 11/05 Started integration of modules
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
library(forecast)
library(verification)
library(doParallel)

ptm <- proc.time() # Start the clock!
cl  <- makeCluster(detectCores())
registerDoParallel(cl)

# From KSSM0100
wm01_00       <- readRDS("0100_import-complete.rds")
importpar     <- readRDS("0100_import-parameter.rds")
s01           <- importpar[1]
s02           <- importpar[2]
s03           <- importpar[3]
sum_of_h      <- importpar[4]
data_size     <- importpar[5]

#===========================================
# Integrated Parameters
#===========================================
cus_list      = seq(1,10)
frontierstp   = 8                       # Number of demand bins
win_size      = c(4,24)                 # Small and large win_size (select only 2)
ahead_t       = seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       = c(0)
in_sample_fr  = 2/3                     # Fraction for diving in- and out-sample
crossvalsize  = 1                       # Number of weeks in the end of in_sample used for crossvalidation
seas_bloc_ws  = 6                       # Number of weeks used for calculating seasonality pattern (6 seems best)
sampling      = 1024                    # For monte-carlo CRPS calculation

#===========================================
# Functions Declarations
#===========================================
fx_evhor <- function (wm01_01,h,in_sample_fr,s02,seas_bloc_ws,crossvalsize){
  event_horizon = ncol(wm01_01)*in_sample_fr+1 + h - s02
  in_sample_ini = event_horizon - min((seas_bloc_ws),(event_horizon %/% s02)) * s02 + 1
  outsample_end = event_horizon + 24
  in_sample_siz = event_horizon - in_sample_ini + 1
  outsample_siz = outsample_end - event_horizon
  inosample_siz = in_sample_siz + outsample_siz
  crossval_ini  = in_sample_siz - crossvalsize * s02 + 1
  event_hrz_mod = crossval_ini - 1
  evhor_out     = c(event_horizon,in_sample_ini,outsample_end,in_sample_siz,outsample_siz,inosample_siz,event_hrz_mod,crossval_ini)
  # starts, in_sample until [7], crossval from [8] to [4], outsample of [5] until [6]
  return(evhor_out)
}

fx_seas   <- function (wm01,s01,s02,sum_of_h,def_evhor){
  wm02    <- foreach (j = 1:nrow(wm01), .combine=c("rbind")) %dopar% {
    wv31i  = wm01[j,1:def_evhor[7]]
    wv32i  = decompose(msts(wv31i,seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
    wv32is = wv32i$seasonal[1:(s02)]
    wv32is
  }
  return(wm02)
}

fx_unseas <- function (wm01,wm02,s02,def_evhor){
  wm04    <- foreach (j = 1:nrow(wm01), .combine=c("rbind")) %dopar% {
    wv33a  = rep(wm02[j,],def_evhor[7]/s02)
    wv33b  = wm02[j,1:(def_evhor[6] - def_evhor[7])]
    wv33   = c(wv33a,wv33b)
    wv34   = wm01[j,] - wv33
    wv34
  }
  return(wm04)
}

fx_fcst_kds <- function (wm04,win_size,def_evhor,sampling){
  densi     <- foreach (j = 1:nrow(wm01)) %dopar% {
    denssmall = density(wm04[j,(def_evhor[7] - win_size[1] + 1):(def_evhor[7])])
    denslarge = density(wm04[j,(def_evhor[7] - win_size[2] + 1):(def_evhor[7])])
    fcstsmall = sample(denssmall$x, sampling, replace=TRUE, prob=denssmall$y)
    fcstlarge = sample(denslarge$x, sampling, replace=TRUE, prob=denslarge$y)
    data.frame(fcstsmall,denssmall$bw,fcstlarge,denslarge$bw)
  }
  return(densi)
}

fx_crps_mc <- function (wm04,densi,def_evhor,sampling){
  crps_mc2 <- foreach (j = 1:nrow(wm01)) %:%
    foreach (i = (def_evhor[4]+1):def_evhor[6], .combine=c("cbind")) %dopar% {
      wv35s = crps(rep(wm04[j,i],sampling),densi[[j]][,1:2])$CRPS
      wv35l = crps(rep(wm04[j,i],sampling),densi[[j]][,3:4])$CRPS
      c(wv35s,wv35l)
    }
  return(crps_mc2)
}






#===========================================
# 
#===========================================
wm01_01    = wm01_00[min(cus_list):length(cus_list),]


def_evhor  = fx_evhor(wm01_01,0,in_sample_fr,s02,seas_bloc_ws,0)
wm01       = wm01_01[,def_evhor[2]:def_evhor[3]]                    # work matrix
wm02       = fx_seas(wm01,s01,s02,sum_of_h,def_evhor)               # in-sample seasonality pattern
wm04       = fx_unseas(wm01,wm02,s02,def_evhor)                     # in-out sample unseasonalised
densi      = fx_fcst_kds(wm04,win_size,def_evhor,sampling)
crps_mc    = fx_crps_mc(wm04,densi,def_evhor,sampling)










#######################################################

#===========================================
# Aggregating demand randomly
#===========================================
cat("\nCreating aggregated demands for random groups of customers\n")
wm01_4D       = list()
wm01_c        = matrix(nrow=(length(cus_list)), ncol=(data_size))
wm01_4Dc <- foreach (c = cus_clu) %dopar% {
  for (i in 1:length(cus_list)) {
    random_pool = round(runif(c,1,length(cus_list)),0)
    cus_samples = wm01_01[random_pool,]
    wm01_c[i,]  = colMeans (cus_samples, na.rm = FALSE, dims = 1)
  }
  wm01_c
}
wm01_4D    = c(list(wm01_01),wm01_4Dc)

#===========================================
# Forecasting and analysis functions
#===========================================
kdscrps <- function(wm01_4D, hrz_lim, win_size, ahead_t, in_sample_fr, crossvalsize, sampling){
  c_group     = seq(1,length(wm01_4D))
  crpsh_CLU   = list()
  wm02_4D   <<- list()
  fcst_4D   <<- list()
  for (c in c_group){
    # cat("\n\n\n----------------------------------------\nGROUP NUMBER ",c," OF ",max(c_group),"\n----------------------------------------")
    wm01    =  wm01_4D[[c]]
    fcsth   =  list()
    crpsh   =  list()
    wm02_3D =  list()
    for (h in hrz_lim) {
      # cat("\nEvent Horizon +",h,"for",nrow(wm01),"customers: ")
      # WorkMatrix 02, 03 & 04 ===================
      # wm02 in-sample seasonality pattern
      # wm03 out-sample load data
      # wm04 in- and out-sample de-seasonalised
      event_horizon = ncol(wm01)*in_sample_fr+1 + h -s02
      in_sample_ini = event_horizon - min(seas_bloc_ws,(event_horizon %/% s02)) * s02 + 1
      outsample_end = (ncol(wm01)*(1-in_sample_fr)) %/% s02 * s02 + event_horizon - (h%/%s02*s02)
      in_sample_siz = event_horizon - in_sample_ini + 1
      outsample_siz = outsample_end - event_horizon
      inosample_siz = in_sample_siz + outsample_siz
      event_hrz_mod = event_horizon - in_sample_ini + 1
      wm02          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=s02)
      wm03          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=outsample_siz)
      wm04          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=(inosample_siz))
      for (j in 1:nrow(wm01_4D[[c]])){
        # if (j == 1){cat("[Seasonality] ")}
        # cat(".")
        wv31i  = wm01[j,(in_sample_ini):(event_horizon)]
        wv31o  = wm01[j,(event_horizon+1):(outsample_end)]
        wv32i  = decompose(msts(wv31i,seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
        wv32is = wv32i$seasonal[1:(s02)]  #####shoud included the /sum_of_h, and 3 lines below?
        wm02[j,1:s02] = wv32is
        wm03[j,1:outsample_siz] = wv31o
        wm04[j,1:inosample_siz] = c(wv31i,wv31o) - rep(wv32is,(inosample_siz)/s02)
      }
      wm02_3D[[match(h,hrz_lim)]] = wm02
      # Forecasting and CRPS =====================
      # cat("[Forecasting]")
      fcstcrpsj <- foreach (j = 1:nrow(wm01_4D[[c]]),
                        .packages=c("forecast","verification")) %dopar% {
        # cat(".")
        crpski = matrix(nrow=length(win_size),ncol=length(ahead_t))
        fcstki = matrix(nrow=length(win_size),ncol=length(ahead_t))
        for (k in win_size){
          for (i in ahead_t){
            densi = density(wm04[j,(event_hrz_mod - k + 1):(event_hrz_mod)])
            fcsti = sample(densi$x, sampling, replace=TRUE, prob=densi$y)
            crpsi = crps(rep(wm04[j,(event_hrz_mod + i)], sampling), data.frame(fcsti, densi$bw))
            crpski[match(k,win_size),i] = crpsi$CRPS
            fcstki[match(k,win_size),i] = mean(fcsti)
          }
        }
        list(fcstki,crpski)
      }
      fcstkij = list()
      crpskij = list()
      for (j in 1:nrow(wm01_4D[[c]])){
        fcstkij[[j]] = fcstcrpsj[[j]][[1]]
        crpskij[[j]] = fcstcrpsj[[j]][[2]]
      }
      crpsh[[match(h,hrz_lim)]]  =  crpskij
      fcsth[[match(h,hrz_lim)]]  =  fcstkij
    }
    fcst_4D[[c]] <<- fcsth
    wm02_4D[[c]] <<- wm02_3D
    # Statistical data from CRPS 4d results =========
    crpsh_All = list()
    crpsh_All_meanki = matrix(nrow=length(win_size),ncol=length(ahead_t))
    crpsh_All_desvki = matrix(nrow=length(win_size),ncol=length(ahead_t))
    # cat('\nCreating statistical data, removing dimention of Horizon and Customers')
    for (k in win_size){
      for (i in ahead_t){
        crpsh_All_tempki = integer(nrow(wm01_4D[[c]])*length(hrz_lim))
        for (h in hrz_lim){
          for (j in 1:nrow(wm01_4D[[c]])){
            crpsh_All_tempki[j+(nrow(wm01_4D[[c]]))*(match(h,hrz_lim)-1)] = crpsh[[match(h,hrz_lim)]][[j]][match(k,win_size),i]
          }
        }
        crpsh_All_meanki[match(k,win_size),i] = mean(crpsh_All_tempki)
        crpsh_All_desvki[match(k,win_size),i] = sd(crpsh_All_tempki)
      }
    }
    crpsh_All = list(crpsh_All_meanki,crpsh_All_desvki)
    crpsh_CLU[[c]] = crpsh_All
  }
  return(crpsh_CLU)
}

ws_crossover_fx <- function(crpsh_CLU){
  ws_crossover  <<- list()
  crpsh_CLU_dyn = crpsh_CLU
  for (j in 1:length(crpsh_CLU)){
    crpsh_CLU_dyn[[j]][[1]] = rbind(crpsh_CLU_dyn[[j]][[1]],0)
    crpsh_CLU_dyn[[j]][[2]] = rbind(crpsh_CLU_dyn[[j]][[2]],0)
    i = 1
    while (crpsh_CLU_dyn[[j]][[1]][1,i] < crpsh_CLU_dyn[[j]][[1]][2,i] &
           i < ncol(crpsh_CLU_dyn[[1]][[1]])){
      crpsh_CLU_dyn[[j]][[1]][3,i] = crpsh_CLU_dyn[[j]][[1]][1,i]
      crpsh_CLU_dyn[[j]][[2]][3,i] = crpsh_CLU_dyn[[j]][[2]][1,i]
      ws_crossover[j] = i+1
      i = i+1
    }
    while (i <= ncol(crpsh_CLU_dyn[[1]][[1]])){
      crpsh_CLU_dyn[[j]][[1]][3,i] = crpsh_CLU_dyn[[j]][[1]][2,i]
      crpsh_CLU_dyn[[j]][[2]][3,i] = crpsh_CLU_dyn[[j]][[2]][2,i]
      i = i+1
    }
  }
  return (crpsh_CLU_dyn)
}

#===========================================
# Function run
#===========================================
crpsh_CLU     = kdscrps(wm01_4D, hrz_lim, win_size, ahead_t, in_sample_fr, crossvalsize, sampling)
crpsh_CLU_dyn = ws_crossover_fx(crpsh_CLU)

print(proc.time() - ptm)        # Stop the clock

#===========================================
# Outputs
#===========================================
saveRDS(crpsh_CLU_dyn,  file="0200_analysis.rds") # crps winsize vs time_adead, mean&desv,                per grouping 
saveRDS(wm01_4D[[1]],   file="0200_origdat1.rds") # data customer vs time_series,                         per grouping
saveRDS(wm02_4D[[1]],   file="0200_seasona1.rds") # seas customer vs seasonal_ts,            per horizon, per grouping
saveRDS(fcst_4D[[1]],   file="0200_forecas1.rds") # fcst winsize vs time_ahead per customer, per horizon, per grouping
saveRDS(parbundl0200,   file="0200_parbundl.rds")
saveRDS(kdscrps,        file="0200_function.rds")
