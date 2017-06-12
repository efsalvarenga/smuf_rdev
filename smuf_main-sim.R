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
# BIG [h] LOOP
#===========================================
bighlpopgr <- list()
bighlpcrps <- list()

for (h in hrz_lim){
  ptm    <- proc.time()
  runkey <- Sys.time()
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  
  #===========================================
  # Individual customers forecast
  #===========================================
  cat("[Ind] ")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06       <- fx_int_fcstgeneric_kdss(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh)
  wv45       <- rowMeans(wl06[[1]])
  sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
  wv46       <- seq(0,frontierstp)^2/frontierstp^2 * sum(wv45)
  
  #===========================================
  # Random groups & evaluation
  #===========================================
  cat("[Rnd] ")
  wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
  wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
  wl06rnd    <- fx_int_fcstgeneric_kdss(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh)
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
                             optgrp   <- genoud(fx_optgrp_sdev, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
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
  res_sdev_kd  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcstgeneric_kdss,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  res_sdev_ag  <- fx_applgrp(optgrp_sdev,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  cat("[OptCVKD] ")
  optgrp_cvkd  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud","foreach"),
                           # .export=c('fcst_mccv'),    #nao precisa mais??
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_crps, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1, print.level=1,
                                                use_arma=F)
                             if(optgrp$value == 10) {
                               grouped = c(rep(0,nrow(wm01_01)))
                             } else {
                               grouped = optgrp$par
                             }
                             grouped
                           }
  bighlpopgr   <- fx_sav_optgrps(c("cvkd",h,frontierstp,length(cus_list),crossvalstps,armalags,crossvalfocus,runkey),optgrp_cvkd)
  res_crps_kd  <- fx_applgrp(optgrp_cvkd,wv46,wm01_01,fx_int_fcstgeneric_kdss,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  res_crps_ag  <- fx_applgrp(optgrp_cvkd,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh,crossvalfocus)
  
  # cat("[OptCVAG] ")
  # optgrp_cvag  <- foreach (i = 1:frontierstp,
  #                          .packages=c("forecast","rgenoud","foreach"),
  #                          .export=c('fcst_mccv'),
  #                          .combine=c("rbind")) %dopar% {
  #                            opt_min_cusd  = wv46[i]
  #                            opt_max_cusd  = wv46[i+1]
  #                            optgrp   <- genoud(fx_optgrp_crps, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
  #                                               starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
  #                                               data.type.int=TRUE,  int.seed=1, print.level=1,
  #                                               use_arma=T)
  #                            if(optgrp$value == 10) {
  #                              grouped = c(rep(0,nrow(wm01_01)))
  #                            } else {
  #                              grouped = optgrp$par
  #                            }
  #                            grouped
  #                          }
  # bighlpopgr   <- fx_sav_optgrps(c("cvag",h,frontierstp,length(cus_list),crossvalstps,armalags,crossvalfocus,runkey),optgrp_cvag)
  # res_crps_ag  <- fx_applgrp(optgrp_cvag,wv46,wm01_01,fx_int_fcstgeneric_armagarch,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)

  bighlpcrps   <- fx_sav_optress(c("sdev+crps_kd+ag",h,frontierstp,length(cus_list),crossvalstps,armalags,crossvalfocus,runkey),
                                 list(c(h,frontierstp,length(cus_list)),cbind(cr01rnd,wv45rnd),res_sdev_kd,res_sdev_ag,res_crps_kd,res_crps_ag))
  # fx_plt_rnd_vs_opt(bighlpcrps[[length(bighlpcrps)]][[2]],c(0.01,0.1),c(0,32),"CRPS")
  cat("\n")
  print(proc.time() - ptm)
}