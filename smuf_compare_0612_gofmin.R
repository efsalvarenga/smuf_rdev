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
cl  <- makeCluster(detectCores())
registerDoParallel(cl)
setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")
savfile = "smuf_compare_0612_gofmin.rds"

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
cus_list      <- seq(1,100)
# frontierstp   <- 5                     # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
win_selec     <- win_size[1]
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(1,9)*667
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(8,8)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
gof.min       <- 0.05                    # GoF crossover value to change ARMA-GARCH to KDS
gof.minseq    <- seq(1,10)^3/10^3

#===========================================
# BIG [h] LOOP Start
#===========================================
resmat.h  <- matrix(nrow=0,ncol=max(ahead_t))
reslis.h  <- rep(list(crpsmat.h),length(gof.minseq)+2)
for (h in hrz_lim){
  ptm <- proc.time()                  # Reset clock
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06kdWS1  <- fx_int_fcstgeneric_kdss(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size[1],is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  wl06kdWS2  <- fx_int_fcstgeneric_kdss(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size[2],is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
  # wl06kd     <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T,armalags,cross_overh)
  wl06ag <- foreach (g = 1:length(gof.minseq)) %do% {
    gof.min = gof.minseq[g]
    fx_int_fcstgeneric_armagarch(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,T,armalags,cross_overh,gof.min)
    }
  for (g in 1:length(gof.minseq)){
    reslis.h[[g]] <- rbind(reslis.h[[g]],colMeans(wl06ag[[g]][[2]],na.rm=T))
  }
  reslis.h[[(length(gof.minseq)+1)]] <-  rbind(reslis.h[[(length(gof.minseq)+1)]],colMeans(wl06kdWS1[[2]],na.rm=T))
  reslis.h[[(length(gof.minseq)+2)]] <-  rbind(reslis.h[[(length(gof.minseq)+2)]],colMeans(wl06kdWS2[[2]],na.rm=T))
  resmat.h  <- matrix(nrow=0,ncol=max(ahead_t))
  for (k in 1:length(reslis.h)){
    resmat.h <- rbind(resmat.h,colMeans(reslis.h[[k]],na.rm=T))
  }
  rownames(resmat.h) <- c(paste("AG",1:length(gof.minseq)),paste("KD",1:2))
  plt.names <- rownames(resmat.h)
  fx_plt_mymat(resmat.h,c(0.08,0.12))
  legend('topright', inset=c(-0.15,-0.2), legend = plt.names,
         lty=1, col=rainbow(length(plt.names)), bty='n', cex=.75, title="Method")
  print(proc.time() - ptm)
  saveRDS(list(resmat.h,reslis.h),  file="smuf_temp_compare.rds")
}
# AG6 was best (gof.min should be ~ 0.2)
#===========================================
# Outputs
#===========================================
saveRDS(list(resmat.h,reslis.h),  file=savfile)
