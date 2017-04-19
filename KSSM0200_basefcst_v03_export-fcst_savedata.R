#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM0200 base forecasting
#-------------------------------------------
# Notes
# 11/04 Modular implemented
# 15/04 Saves forecasts and data for later use
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

# source('./KSSM0100_import_v01.R')

#===========================================
# Parameters & Start data
#===========================================
# From KSSM0100
sum_of_h      = sum_of_h
data_size     = data_size

# New parameters
cus_list      = seq(1,200)
cus_clu       = c(10,50,100,200)
wm01_01       = wm01_00[min(cus_list):length(cus_list),]
win_size      = c(4,24)
ahead_t       = seq(1, (12/sum_of_h))
hrz_lim       = c(0,846,1692,2538) #seq(0,(150/sum_of_h)) * 23
in_sample_fr  = 2/3           # Fraction for diving in- and out-sample
sampling      = 1024          # For monte-carlo CRPS calculation

# Parameter Bundle
periodpar     = c(s01,s02,s03,sum_of_h,data_size,in_sample_fr)
parbundl0200  = list(periodpar,cus_list,cus_clu,win_size,ahead_t,hrz_lim,sampling)

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
    wm01_c[i,] = colMeans (cus_samples, na.rm = FALSE, dims = 1)
  }
  wm01_c
}
wm01_4D    = c(list(wm01_01),wm01_4Dc)

#===========================================
# Forecasting and analysis function
#===========================================
kdscrps <- function(wm01_4D, hrz_lim, win_size, ahead_t, sampling){
  c_group     = seq(1,length(wm01_4D))
  crpsh_CLU   = list()
  wm02_4D   <<- list()
  fcst_4D   <<- list()
  for (c in c_group){
    cat("\n\n\n----------------------------------------\nGROUP NUMBER ",c," OF ",max(c_group),"\n----------------------------------------")
    wm01    =  wm01_4D[[c]]
    fcsth   =  list()
    crpsh   =  list()
    wm02_3D =  list()
    for (h in hrz_lim) {
      cat("\nEvent Horizon +",h,"for",length(cus_list),"customers: ")
      # WorkMatrix 02, 03 & 04 ===================
      # wm02 in-sample seasonality pattern
      # wm03 out-sample load data
      # wm04 in- and out-sample de-seasonalised
      event_horizon = data_size*in_sample_fr+1 + h
      in_sample_ini = event_horizon - event_horizon %/% s02 * s02 + 1
      outsample_end = (data_size*(1-in_sample_fr)) %/% s02 * s02 + event_horizon - s02 - (h%/%s02*s02)
      in_sample_siz = event_horizon - in_sample_ini + 1
      outsample_siz = outsample_end - event_horizon
      inosample_siz = in_sample_siz + outsample_siz
      event_hrz_mod = event_horizon - in_sample_ini + 1
      wm02          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=s02)
      wm03          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=outsample_siz)
      wm04          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=(inosample_siz))
      for (j in 1:nrow(wm01_4D[[c]])){
        if (j == 1){cat("[Seasonality] ")}
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
      cat("[Forecasting]")
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
      for (j in 1:length(cus_list)){
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
    cat('\nCreating statistical data, removing dimention of Horizon and Customers')
    for (k in win_size){
      for (i in ahead_t){
        crpsh_All_tempki = integer(length(cus_list)*length(hrz_lim))
        for (h in hrz_lim){
          for (j in 1:nrow(wm01_4D[[c]])){
            crpsh_All_tempki[j+(length(cus_list))*(match(h,hrz_lim)-1)] = crpsh[[match(h,hrz_lim)]][[j]][match(k,win_size),i]
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

#===========================================
# Function run
#===========================================
crpsh_CLU = kdscrps(wm01_4D, hrz_lim, win_size, ahead_t, sampling)

ws_crossover  = list()
crpsh_CLU_dyn = crpsh_CLU
for (j in 1:length(crpsh_CLU)){
  crpsh_CLU_dyn[[j]][[1]] = rbind(crpsh_CLU_dyn[[j]][[1]],0)
  crpsh_CLU_dyn[[j]][[2]] = rbind(crpsh_CLU_dyn[[j]][[2]],0)
  i = 1
  while (crpsh_CLU_dyn[[j]][[1]][1,i] < crpsh_CLU_dyn[[j]][[1]][2,i]){
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

print(proc.time() - ptm)        # Stop the clock

#===========================================
# Outputs
#===========================================
saveRDS(crpsh_CLU_dyn,  file="0200_analysis.rds")
saveRDS(wm01_4D[[1]],   file="0200_origdat1.rds")
saveRDS(wm02_4D[[1]],   file="0200_seasona1.rds")
saveRDS(fcst_4D[[1]],   file="0200_forecas1.rds")
saveRDS(parbundl0200,   file="0200_parbundl.rds")
saveRDS(kdscrps,        file="0200_function.rds")