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
savfile = "smuf_run_0720_defheur_clulight_int02.rds"

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
cus_list      <- seq(1,200)
frontierstp   <- 15             # Number of demand bins (Stepwise frontier for portfolio optimisation)
frontierexp   <- 1                     # Exponentiality of frontier steps
max.gen       <- 300                     # For genetic opt
waitgen       <- 50                      # For genetic opt
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
win_selec     <- win_size[2]
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1,24)               # Up to s02
hrz_lim       <- seq(0,10)*113*3            # Rolling forecasts steps {seq(0:167)*113} is comprehensive
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 8                      # Steps used for multiple crossvalidation (Only KDE)
crossvalfocus <- c(1)                  # What period is focused when running crossvalidation
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(5,5)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
gof.min       <- 0.2                    # GoF crossover value to change ARMA-GARCH to KDS

#===========================================
# Call simulator
#===========================================
bighlpopgr <- list()
bighlpcrps <- list()
myleg = c("Random","sdev","seas non-scaled all","seas non-scaled all02","seas non-scaled all05","seas non-scaled all08")

for (h in hrz_lim){
  ptm    <- proc.time()
  runkey <- Sys.time()
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  
  #===========================================
  # Individual customers forecast
  #===========================================
  cat("[Ind] ")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06       <- fx_int_fcstgeneric_kdss(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh,gof.min)
  wv45       <- rowMeans(wl06[[1]])
  sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
  wv46       <- seq(0,frontierstp)^frontierexp/frontierstp^frontierexp * sum(wv45)
  
  #===========================================
  # Random groups & evaluation
  #===========================================
  cat("[Rnd] ")
  wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
  wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
  wl06rnd    <- fx_int_fcstgeneric_kdss(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh,gof.min)
  wv45rnd    <- as.numeric(rowMeans(wl06rnd[[1]]) * rowSums(wm01_02l[[2]]))
  # sd01rnd    <- as.numeric(fx_sd_mymat(wl06rnd[[3]]))
  cr01rnd    <- rowMeans(wl06rnd[[2]])
  
  #===========================================
  # Optimising groups & evaluation
  #===========================================
  cat("[OptSDEV] ")
  optgrp_sdev  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_sdev, nvars=nrow(wm01_01), max.generations=max.gen, wait.generations=waitgen,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  bighlpopgr   <- fx_sav_optgrps(c("sdev",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_sdev)
  res_sdev_ag  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  wm13n        <- t(scale(t(wl06[[4]]))) # for clulight03, feed optSDEV as initial
  
  cat("[OptSEAF_na] ")
  wm13seaf     <- as.matrix(wl06[[4]])
  optgrp_sfna  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_seaf, nvars=nrow(wm01_01), max.generations=max.gen, wait.generations=waitgen,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  bighlpopgr   <- fx_sav_optgrps(c("sfna",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_sfna)
  res_sfna_ag  <- fx_applgrp(optgrp_sfna,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  cat("[OptSEAF_na02] ")
  wm13seaf     <- as.matrix(wl06[[4]])
  optgrp_sfna02  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_ssmix, nvars=nrow(wm01_01), max.generations=max.gen, wait.generations=waitgen,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1,
                                                defratsd=0.2)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  bighlpopgr   <- fx_sav_optgrps(c("sfna02",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_sfna02)
  res_sfna02_ag  <- fx_applgrp(optgrp_sfna02,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  cat("[OptSEAF_na05] ")
  wm13seaf     <- as.matrix(wl06[[4]])
  optgrp_sfna05  <- foreach (i = 1:frontierstp,
                             .packages=c("forecast","rgenoud"),
                             .combine=c("rbind")) %dopar% {
                               opt_min_cusd  = wv46[i]
                               opt_max_cusd  = wv46[i+1]
                               optgrp   <- genoud(fx_optgrp_ssmix, nvars=nrow(wm01_01), max.generations=max.gen, wait.generations=waitgen,
                                                  starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                  data.type.int=TRUE,  int.seed=1,
                                                  print.level=1,
                                                  defratsd=0.5)
                               if(optgrp$value == 10) {
                                 grouped = c(rep(0,nrow(wm01_01)))
                               } else {
                                 grouped = optgrp$par
                               }
                               grouped
                             }
  bighlpopgr   <- fx_sav_optgrps(c("sfna05",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_sfna05)
  res_sfna05_ag  <- fx_applgrp(optgrp_sfna05,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  cat("[OptSEAF_na08] ")
  wm13seaf     <- as.matrix(wl06[[4]])
  optgrp_sfna08  <- foreach (i = 1:frontierstp,
                             .packages=c("forecast","rgenoud"),
                             .combine=c("rbind")) %dopar% {
                               opt_min_cusd  = wv46[i]
                               opt_max_cusd  = wv46[i+1]
                               optgrp   <- genoud(fx_optgrp_ssmix, nvars=nrow(wm01_01), max.generations=max.gen, wait.generations=waitgen,
                                                  starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                  data.type.int=TRUE,  int.seed=1,
                                                  print.level=1,
                                                  defratsd=0.8)
                               if(optgrp$value == 10) {
                                 grouped = c(rep(0,nrow(wm01_01)))
                               } else {
                                 grouped = optgrp$par
                               }
                               grouped
                             }
  bighlpopgr   <- fx_sav_optgrps(c("sfna08",h,frontierstp,length(cus_list),crossvalstps,armalags,runkey),optgrp_sfna08)
  res_sfna08_ag  <- fx_applgrp(optgrp_sfna08,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  bighlpcrps   <- fx_sav_optress(c("sdev+seafs_mix",h,frontierstp,length(cus_list),crossvalstps,armalags,crossvalfocus,runkey),
                                   list(c(h,frontierstp,length(cus_list)),cbind(cr01rnd,wv45rnd),res_sdev_ag,res_sfna_ag,res_sfna02_ag,res_sfna05_ag,res_sfna08_ag))
  
  saveRDS(list(bighlpopgr,bighlpcrps),  file=savfile)
  
  fx_plt_rnd_vs_opt(bighlpcrps[[length(bighlpcrps)]][[2]],c(0,0.1),c(0,sum(wv45)),myleg,"CRPS")
  cat("\n")
  print(proc.time() - ptm)
}

biglpcrpsavg = list()

for (j in 1:7){
  biglpcrpsavg[[j]]=bighlpcrps[[1]][[2]][[j]]
}
for (i in 2:length(bighlpcrps)){
  for (j in 1:7){
    biglpcrpsavg[[j]] = biglpcrpsavg[[j]] + bighlpcrps[[i]][[2]][[j]]
  }
}
for (j in 1:7){
  biglpcrpsavg[[j]]=biglpcrpsavg[[j]]/11
}
fx_plt_rnd_vs_opt(biglpcrpsavg,c(0,0.1),c(0,sum(wv45)),myleg,"CRPS")

for (j in 3:7){
  cat('\n',myleg[(j-1)],' ',mean(biglpcrpsavg[[j]][,1],na.rm=T))
}

saveRDS(biglpcrpsavg,  file=paste("summary_",savfile,sep=""))
