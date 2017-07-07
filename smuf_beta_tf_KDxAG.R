# TF to decide KD or AG

#===========================================
# Initialising
#===========================================
setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")
cl  <- makeCluster(detectCores())
registerDoParallel(cl)
savfile = "smuf_compare_0629_KDxAG_large.rds"

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
cus_list      <- seq(1,10)
# frontierstp   <- 5                     # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
win_selec     <- win_size[2]
# cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
# ahead_t       <- seq(1, (60/sum_of_h))   # Up to s02
h             <- 0
# hrz_lim       <- seq(20,50)*113            # Rolling forecasts steps {seq(0:167)*113} is comprehensive
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
# crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
# crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(5,5)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
gof.min       <- 0                    # GoF crossover value to change ARMA-GARCH to KDS




#===========================================
# 
#===========================================
wm01_01    <- wm01_00[min(cus_list):max(cus_list),]

out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,0)
wm01       <- wm01_01[,out_evhor[2]:out_evhor[3]]                         # work matrix
wl02       <- fx_seas2(wm01,s01,s02,sum_of_h,out_evhor)                   # in-sample seasonality pattern (s,r,t)
wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]                        # out-sample original load data
wm13       <- fx_unseas2(wm01,wl02,s02,out_evhor)                         # out-sample estimated trend + seas
wm14       <- wl02[[2]]                                                   # in-sample noise




fcst_mc    <- fx_fcst_armagarchs(wm14,armalags,win_selec,ahead_t,out_evhor,sampling,cross_overh,gof.min) # returns list with next ahead_t fcst and sd
wm03fcst   <- fx_fcstgeneric(fcst_mc,out_evhor,wm13)
wm05       <- fx_crpsgeneric(wm03,wm13,wm14,fcst_mc,out_evhor,sampling)

