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
savfile = "smuf_run_0611_1400.rds"

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
frontierstp   <- 10                      # Number of demand bins (Stepwise frontier for portfolio optimisation)
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
# Call simulator
#===========================================s
source("smuf_main.R")
saveRDS(list(bighlpopgr,bighlpcrps),  file=savfile)