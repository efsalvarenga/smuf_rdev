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
source("smuf_main-fxs.R")
cl  <- makeCluster(detectCores())
registerDoParallel(cl)
savfile = "smuf_compare_0831_KD_mult-win_kdss1.rds"

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
# change for superun: cus 1000, hrz_lim steps, gof steps, arma-lags 10
cus_list      <- seq(1,50)#300)
# frontierstp   <- 5                     # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
win_selec     <- win_size[2]
# cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (72/sum_of_h))   # Up to s02
hrz_lim       <- seq(1:180)*37            # Rolling forecasts steps {seq(0:167)*113} is comprehensive
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
# crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
# crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(5,5)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
gof.min       <- 0.1                    # GoF crossover value to change ARMA-GARCH to KDS

#===========================================
# BIG [h] LOOP Start
#===========================================
resmat.h  <- matrix(nrow=0,ncol=max(ahead_t))
reslis.h  <- rep(list(resmat.h),5)
wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
for (h in hrz_lim){
  ptm <- proc.time()                  # Reset clock
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  cat("KD04 ")
  wl06kdWS04  <- fx_int_fcstgeneric_kdss2(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,4,is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  cat("KD06 ")
  wl06kdWS06  <- fx_int_fcstgeneric_kdss2(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,6,is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  cat("KD12 ")
  wl06kdWS12  <- fx_int_fcstgeneric_kdss2(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,12,is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  cat("KD24 ")
  wl06kdWS24  <- fx_int_fcstgeneric_kdss2(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,24,is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  cat("KD48 ")
  wl06kdWS48  <- fx_int_fcstgeneric_kdss2(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,48,is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  cat(" [Consolidating]\n")
  reslis.h[[1]] <- rbind(reslis.h[[1]],colMeans(wl06kdWS04[[2]],na.rm=T))
  reslis.h[[2]] <- rbind(reslis.h[[2]],colMeans(wl06kdWS06[[2]],na.rm=T))
  reslis.h[[3]] <- rbind(reslis.h[[3]],colMeans(wl06kdWS12[[2]],na.rm=T))
  reslis.h[[4]] <- rbind(reslis.h[[4]],colMeans(wl06kdWS24[[2]],na.rm=T))
  reslis.h[[5]] <- rbind(reslis.h[[5]],colMeans(wl06kdWS48[[2]],na.rm=T))
  resmat.h  <- matrix(nrow=0,ncol=max(ahead_t))
  for (k in 1:length(reslis.h)){
    resmat.h <- rbind(resmat.h,colMeans(reslis.h[[k]],na.rm=T))
  }
  rownames(resmat.h) <- c("KD04","KD06","KD12","KD24","KD48")
  plt.names <- rownames(resmat.h)
  fx_plt_mymat(resmat.h,c(0.06,0.13))
  legend('topright', inset=c(0,0), legend = plt.names,
         lty=1, col=rainbow(length(plt.names)), bty='n', cex=.75, title="Method")
  print(proc.time() - ptm)
  saveRDS(list(resmat.h,reslis.h),  file="smuf_temp_compare.rds")
}
#===========================================
# Outputs
#===========================================
saveRDS(list(resmat.h,reslis.h),  file=savfile)
