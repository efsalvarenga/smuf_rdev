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
source("smuf_main-fxs.R")
savfile = "smuf_run_0704_sd04_med.rds"

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
#cus_list to 1000, stp to 150 (detectcores), hrz_lim larger (0:167)*113), turn on CV
cus_list      <- seq(1,100)
frontierstp   <- 16                      # Number of demand bins (Stepwise frontier for portfolio optimisation)
frontierexp   <- 1.2                     # Exponentiality of frontier steps
max.gen       <- 300                     # For genetic opt
waitgen       <- 50                      # For genetic opt
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
win_selec     <- win_size[2]
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (12/sum_of_h))    # Up to s02
hrz_lim       <- 0 #seq(5,6)*113            # Rolling forecasts steps {seq(0:167)*113} is comprehensive
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 16                      # Steps used for multiple crossvalidation (Only KDE)
crossvalfocus <- c(4)                  # What period is focused when running crossvalidation
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(5,5)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
gof.min       <- 0.2                    # GoF crossover value to change ARMA-GARCH to KDS

#===========================================
# Call simulator
#===========================================
OptCVKD = F
OptCVAG = F
source("smuf_main-optgrp.R")
saveRDS(list(bighlpopgr,bighlpcrps),  file=savfile)
