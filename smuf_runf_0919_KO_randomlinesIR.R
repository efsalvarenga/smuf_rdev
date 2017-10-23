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
  # library(tidyr)
  # library(dplyr)
  source("smuf_main-fxs.R")
  savfile = "smuf_runf_0919_KO_randomlinesIR.rds"
  
  wm01_00       <- readRDS("smuf_import-completeIRhour.rds")
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
  frontierstp   <- 16             # Number of demand bins (Stepwise frontier for portfolio optimisation)
  frontierexp   <- 1                     # Exponentiality of frontier steps
  max.gen       <- 100                     # For genetic opt
  waitgen       <- 10                      # For genetic opt
  win_size      <- c(4,24)                 # Small and large win_size (select only 2)
  win_selec     <- win_size[2]
  cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
  ahead_t       <- seq(1,72)               # Up to s02
  hrz_lim       <- seq(1,167)*29            # Rolling forecasts steps {seq(0:167)*113} is comprehensive
  in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
  crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
  crossvalstps  <- 16                      # Steps used for multiple crossvalidation (Only KDE)
  crossvalfocus <- c(1)                  # What period is focused when running crossvalidation
  is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
  sampling      <- 1024                    # For monte-carlo CRPS calculation
  armalags      <- c(5,5)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
  gof.min       <- 0.05                    # GoF crossover value to change ARMA-GARCH to KDS
  
  #===========================================
  # Call simulator
  #===========================================
  bigrndno <- data.frame(V1=double(),
                        hor=character(),
                        CRPS=double())
  rndres_big <- list()
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
    # wv45       <- rowMeans(wl06[[1]])
    # sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
    # wv46       <- seq(0,frontierstp)^frontierexp/frontierstp^frontierexp * sum(wv45)
    wv51         <- colMeans(wl06[[2]])
    rndres       <- wv51
    #===========================================
    # Random groups & evaluation
    #===========================================
    cat("[Rnd] ")
    rnd.names <- c(2,3,4,5,10,20,30,40,50,100,150,200)
    for (c in rnd.names){
      cat(c," ")
      matmult <- matrix(0,length(cus_list),length(cus_list))
      for (j in 1:length(cus_list)){
        vecmult <- rep(0,length(cus_list))
        vecmult[sample(cus_list,c,replace=F)] = 1/c
        matmult[j,] = vecmult
      }
      wm01_02    <- matmult %*% wm01_01
      wl06_02    <- fx_int_fcstgeneric_kdss(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh,gof.min)
      rndres     <- rbind(rndres,colMeans(wl06_02[[2]]))
    }
    # # Previous Rndgrp implementation [slow]
    # wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
    # wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
    # wl06rnd    <- fx_int_fcstgeneric_kdss(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_selec,is_wins_weeks,crossvalsize,fcst_run,armalags,cross_overh,gof.min)
    # wv45rnd    <- as.numeric(rowMeans(wl06rnd[[1]]) * rowSums(wm01_02l[[2]]))
    # sd01rnd    <- as.numeric(fx_sd_mymat(wl06rnd[[3]]))
    # cr01rnd    <- rowMeans(wl06rnd[[2]])
    # rnd_per_no <- as.data.frame(cbind(as.numeric(rowSums(wm01_02l[[2]])),wl06rnd[[2]]))
    # rndnosimpl <- rnd_per_no %>%
    #   gather(hor,crps,-V1) %>%
    #   group_by(V1,hor) %>%
    #   summarise(CRPS=mean(crps))
    # bigrndno   <- rbind(bigrndno,as.data.frame(rndnosimpl))
    # saveRDS(bigrndno,  file=savfile)
    # 
    # plt.names <- c(1,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200)
    # mymat <- as.matrix(
    #   bigrndno %>%
    #     group_by(V1,hor) %>%
    #     summarise(CRPS=mean(CRPS)) %>%
    #     filter(V1 %in% plt.names) %>%
    #     spread(hor,CRPS)
    #   )
    fx_plt_mymat(rndres,c(0,0.2))
    legend('topright', inset=c(0,0), legend = c(1,rnd.names),
           lty=1, col=rainbow(1+length(rnd.names)), bty='n', cex=.75, title="Method")    
    rndres_big[[match(h,hrz_lim)]] <- rndres
    cat("\n")
    print(proc.time() - ptm)
  }
  rndres_sum <- Reduce("+", rndres_big) / length(rndres_big)
  fx_plt_mymat(rndres_sum,c(0,0.5))
  legend('topright', inset=c(0,0), legend = c(1,rnd.names),
         lty=1, col=rainbow(1+length(rnd.names)), bty='n', cex=.75, title="Method")    
  
  saveRDS(rndres_sum,'rndres_sum.rds')
  