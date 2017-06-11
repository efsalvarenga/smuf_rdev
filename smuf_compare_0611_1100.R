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
# Initialising
#===========================================
setwd("~/GitRepos/smuf_rdev")
source("smuf_fxs.R")
savfile = "smuf_compare_0611_1100.rds"

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
# frontierstp   <- 5                     # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(1,12)*667
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(8,8)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)

#===========================================
# BIG [h] LOOP Start
#===========================================
crpskdmath <- matrix(nrow=0,ncol=max(ahead_t))
crpsagmath <- matrix(nrow=0,ncol=max(ahead_t))
crpskdmatc <- matrix(nrow=0,ncol=(max(cus_list)-min(cus_list)+1))
crpsagmatc <- matrix(nrow=0,ncol=(max(cus_list)-min(cus_list)+1))
for (h in hrz_lim){
  ptm <- proc.time()                  # Reset clock
  cl  <- makeCluster(detectCores())   # Reset workers
  registerDoParallel(cl)
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06kd     <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T,armalags,cross_overh)
  wl06ag     <- fx_int_fcstgeneric_armagarch(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh)
  crpskdmath <- rbind(crpskdmath,colMeans(wl06kd[[2]],na.rm=T))
  crpsagmath <- rbind(crpsagmath,colMeans(wl06ag[[2]],na.rm=T))
  crpskdmatc <- rbind(crpskdmatc,rowMeans(wl06kd[[2]],na.rm=T))
  crpsagmatc <- rbind(crpsagmatc,rowMeans(wl06ag[[2]],na.rm=T))
  cat("Average CRPS for KDS:",mean(wl06kd[[2]])," for ARMA-GARCH:",mean(wl06ag[[2]],na.rm=T),"\n")
  fx_plt_mymat(rbind(colMeans(crpskdmath),colMeans(crpsagmath)),c(0.05,0.15))
  legend('topright', inset=c(0,0), legend = c("KDS","ARMA-GARCH"),
         lty=1, col=rainbow(2), bty='n', cex=.75, title="Grouping")
  fx_plt_mymat(rbind(colMeans(crpskdmatc),colMeans(crpsagmatc)),c(0.0,0.2))
  legend('topright', inset=c(0,0), legend = c("KDS","ARMA-GARCH"),
         lty=1, col=rainbow(2), bty='n', cex=.75, title="Grouping")
  print(proc.time() - ptm)
  crpscompare = list(crpskdmath,crpsagmath,crpskdmatc,crpsagmatc)
  saveRDS(crpscompare,  file="smuf_compare-temp.rds")
}

#===========================================
# Outputs
#===========================================
saveRDS(crpscompare,  file=savfile)
