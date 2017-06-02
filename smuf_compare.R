#===========================================
# Smart Metering Uncertainty Forecasting
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# smuf_compare arma-garch vs kds [individual]
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
library(forecast)
library(verification)
library(doParallel)
library(rgenoud)
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
cus_list      <- seq(1,50)
# frontierstp   <- 5                       # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(0,2)*2069
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
seas_bloc_ws  <- 6                       # Number of weeks used for calculating seasonality pattern (6 seems best)
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(6,6)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)

#===========================================
# Functions Declarations
#===========================================

#===========================================
# BIG [h] LOOP Start
#===========================================
crpskdmath <- matrix(nrow=0,ncol=max(ahead_t))
crpsagmath <- matrix(nrow=0,ncol=max(ahead_t))
crpskdmatc <- matrix(nrow=0,ncol=max(cus_list))
crpsagmatc <- matrix(nrow=0,ncol=max(cus_list))
for (h in hrz_lim){
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  wm01_01    <- wm01_00[min(cus_list):length(cus_list),]
  wl06kd     <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,seas_bloc_ws,crossvalsize,T)
  wl06ag     <- fx_int_fcstgeneric_armagarch(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,seas_bloc_ws,crossvalsize)
  crpskdmath <- rbind(crpskdmath,colMeans(wl06kd[[2]]))
  crpsagmath <- rbind(crpsagmath,colMeans(wl06ag[[2]]))
  crpskdmatc <- rbind(crpskdmatc,rowMeans(wl06kd[[2]]))
  crpsagmatc <- rbind(crpsagmatc,rowMeans(wl06ag[[2]]))
  print(proc.time() - ptm)
}

# fx_plt_mymat(crpskdmath,c(0.05,0.5))
# fx_plt_mymat(crpsagmath,c(0.05,0.5))
# fx_plt_mymat(crpskdmatc,c(0.05,0.5))
# fx_plt_mymat(crpsagmatc,c(0.05,0.5))

crpscompare = list(crpskdmath,crpsagmath,crpskdmatc,crpsagmatc)

#===========================================
# Outputs
#===========================================
saveRDS(crpscompare,  file="smuf_compare-kd_ag.rds")
