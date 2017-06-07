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
cus_list      <- seq(1,30)
frontierstp   <- 5                       # Number of demand bins (Stepwise frontier for portfolio optimisation)
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(0,1)*2069
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
crossvalstps  <- 2                       # Steps used for multiple crossvalidation (Only KDE)
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(2,3)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)

#===========================================
# Functions Declarations: Modules
#===========================================
fx_evhor <- function (wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,crossvalsize){
  event_horizon <- ncol(wm01_01) * in_sample_fr + 1 + h - s02
  in_sample_ini <- event_horizon - min((is_wins_weeks),(event_horizon %/% s02)) * s02 + 1
  outsample_end <- event_horizon + max(ahead_t)
  in_sample_siz <- event_horizon - in_sample_ini + 1
  outsample_siz <- outsample_end - event_horizon
  inosample_siz <- in_sample_siz + outsample_siz
  crossval_ini  <- in_sample_siz - crossvalsize * s02 + 1
  event_hrz_mod <- crossval_ini - 1
  if (crossvalsize==0){
    crossval_ini=NA
    crossval_end=NA
  }
  evhor_out     <- c(event_horizon,in_sample_ini,outsample_end,in_sample_siz,outsample_siz,inosample_siz,event_hrz_mod,crossval_ini)
  return(evhor_out)
}

fx_seas   <- function (wm01,s01,s02,sum_of_h,def_evhor){
  wm02    <- foreach (j = 1:nrow(wm01), .combine=c("rbind"), .packages=c("forecast")) %dopar% {
    wv31i  <- wm01[j,1:def_evhor[7]]
    wv32i  <- decompose(msts(wv31i,seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
    wv32is <- wv32i$seasonal[1:(s02)]
    wv32is
  }
  return(wm02)
}

fx_unseas <- function (wm01,wm02,s02,def_evhor){
  wm04    <- foreach (j = 1:nrow(wm01), .combine=c("rbind")) %dopar% {
    wv33a  <- rep(wm02[j,],def_evhor[7]/s02)
    wv33b  <- wm02[j,1:(def_evhor[6] - def_evhor[7])]
    wv33   <- c(wv33a,wv33b)
    wv34   <- wm01[j,] - wv33
    wv34
  }
  return(wm04)
}

fx_seas2  <- function (wm01,s01,s02,sum_of_h,def_evhor){
  wm12    <- foreach (j = 1:nrow(wm01), .combine=c("rbind"), .packages=c("forecast")) %dopar% {
    wv31i  <- wm01[j,1:def_evhor[7]]
    wv32i  <- stl(msts(wv31i,seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)), s.window="periodic", robust=TRUE)
    wv32is <- wv32i$time.series[1:s02,1]
    wv32it <- wv32i$time.series[,2]
    wv32ir <- wv32i$time.series[,3]
    c(wv32is,wv32ir,wv32it)
  }
  sepp  <- s02+(ncol(wm12)-s02)/2
  wm12s <- wm12[,1:s02]
  wm12r <- wm12[,(s02+1):sepp]
  wm12t <- wm12[,(sepp+1):ncol(wm12)]
  return(list(wm12s,wm12r,wm12t))
}

fx_unseas2 <- function (wm01,wm02,s02,def_evhor){
  wm14     <- foreach (j = 1:nrow(wm01), .combine=c("rbind")) %dopar% {
    wv33b  <- wm02[[1]][j,1:(def_evhor[6] - def_evhor[7])]
    wv33c  <- rep(mean(tail(wm02[[3]][j,],1)),(def_evhor[6] - def_evhor[7]))
    wv33b + wv33c
  }
  return(wm14)
}

fx_fcst_kds <- function (wm04,win_size,def_evhor,sampling){
  fcst_mc2    <- foreach (j = 1:nrow(wm04)) %dopar% {
    denssmall <- density(wm04[j,(def_evhor[7] - win_size[1] + 1):(def_evhor[7])])
    denslarge <- density(wm04[j,(def_evhor[7] - win_size[2] + 1):(def_evhor[7])])
    fcstsmall <- sample(denssmall$x, sampling, replace=TRUE, prob=denssmall$y)
    fcstlarge <- sample(denslarge$x, sampling, replace=TRUE, prob=denslarge$y)
    data.frame(fcstsmall,denssmall$bw,fcstlarge,denslarge$bw)
  }
  return(fcst_mc2)
}

fx_fcst_kds_quickvector <- function (runvec,win_size,def_evhor,sampling,cross_overh){
  denssmall  <- density(runvec[(def_evhor[7] - win_size[1] + 1):(def_evhor[7])])
  denslarge  <- density(runvec[(def_evhor[7] - win_size[2] + 1):(def_evhor[7])])
  fcstsmall  <- sample(denssmall$x, sampling, replace=TRUE, prob=denssmall$y)
  fcstlarge  <- sample(denslarge$x, sampling, replace=TRUE, prob=denslarge$y)
  kdsim.mean <- rbind(t(replicate(cross_overh,fcstsmall)),t(replicate((max(ahead_t)-cross_overh),fcstlarge)))
  kdsim.desv <- rbind(t(replicate(cross_overh,rep(denssmall$bw,sampling))),t(replicate((max(ahead_t)-cross_overh),rep(denslarge$bw,sampling))))
  fcst_kds_qv <- list(kdsim.mean,kdsim.desv)
  return(fcst_kds_qv)
}

fx_fcst_armagarch <- function (wm14,armalags,win_size,ahead_t,out_evhor,sampling,cross_overh){
  fcst_armagarch <- foreach (j = 1:nrow(wm14), .packages=c("rugarch"), .export='fx_fcst_kds_quickvector') %dopar% {
    runvec       <- wm14[j,]
    # Defining ARMA lags
    final.bic <- matrix(nrow=0,ncol=4)
    for (p in 0:armalags[1]) for (q in 0:armalags[2]) {
      if ( p == 0 && q == 0) {
        next
      }
      arimaFit = tryCatch(arima(runvec, order=c(p, 0, q)),
                          error=function(err) FALSE,
                          warning=function( err ) FALSE )
      if( !is.logical(arimaFit) ) {
        final.bic <- rbind(final.bic,c(p,q,AIC(arimaFit,k=log(out_evhor[7])),AIC(arimaFit)))
      } else {
        next
      }
    }
    final.ord <- final.bic[sort.list(final.bic[,3]), ]
    if (nrow(final.ord)==0) {             # if no ARMA(p,q) fits, go with kde_quickvector
      simdata <- fx_fcst_kds_quickvector(runvec,win_size,out_evhor,sampling,cross_overh)
    } else {                              # fitting ARMA-GARCH parameters and simulating
      fit        = F
      final.ordl = 0
      while (!is.logical(fit) == F) {
        spec = ugarchspec(variance.model=list(garchOrder=c(1,1)),
                          mean.model=list(armaOrder=c(final.ord[(final.ordl+1),1], final.ord[(final.ordl+1),2]), include.mean=T),
                          distribution.model="sged")
        fit = tryCatch(ugarchfit(spec, runvec, solver = 'hybrid'),
                       error=function(err) FALSE, warning=function(err) FALSE)
        if(final.ordl <= nrow(final.ord)) {final.ordl = final.ordl+1}
        next
      }
      if (gof(fit, groups=c(2016))[3] >= 0.05) {
        sim1    <- ugarchsim(fit, n.sim = max(ahead_t), m.sim = sampling)
        simdata <- list(sim1@simulation$seriesSim,sim1@simulation$sigmaSim)
      } else {
        simdata <- fx_fcst_kds_quickvector(runvec,win_size,out_evhor,sampling,cross_overh)
      }
    }
    simdata
  }
  return(fcst_armagarch)
}

fx_crps_mc <- function (wm04,fcst_mc,def_evhor,sampling){
  crps_mc2 <- foreach (j = 1:nrow(wm04)) %:%
    foreach (i = (def_evhor[4]+1):def_evhor[6], .combine=c("cbind"), .packages=c("verification")) %dopar% {
      wv35s <- crps(rep(wm04[j,i],sampling),fcst_mc[[j]][,1:2])$CRPS
      wv35l <- crps(rep(wm04[j,i],sampling),fcst_mc[[j]][,3:4])$CRPS
      c(wv35s,wv35l)
    }
  return(crps_mc2)
}

fx_crpsgeneric <- function (wm03,wm13,wm14,fcst_mc,def_evhor,sampling){
  crps_mc2 <- foreach (j = 1:nrow(wm14), .combine=c("rbind")) %:%
    foreach (i = 1:ncol(wm03), .combine=c("cbind"), .packages=c("verification")) %dopar% {
      wv35 <- crps(rep((wm03[j,i]-wm13[j,i]),sampling),data.frame(fcst_mc[[j]][[1]][i,],fcst_mc[[j]][[2]][i,]))$CRPS
      c(wv35)
    }
  return(crps_mc2)
}

fx_crossover <- function(fcst_mc,crps_mc,wm02,def_evhor){
  wscoj <- foreach (j = 1:nrow(wm02),.combine=c("cbind")) %dopar% {
    i=1
    wscrossover <- 1
    while (match(min(crps_mc[[j]][,i]),crps_mc[[j]][,i]) < 2 & i < def_evhor[5]){
      wscrossover <- i+1
      i <- i+1
    }
    wscrossover
  }
  return(wscoj)
}

fx_getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

fx_fcst_wm <- function(fcst_mc,cvcojmean,out_evhor,wm02){
  fcst_co <- foreach (j = 1:nrow(wm02),.combine=c("rbind")) %dopar% {
    wv36f <- matrix(ncol=out_evhor[5])
    fcst1 <- colMeans(fcst_mc[[j]])[1]
    fcst2 <- colMeans(fcst_mc[[j]])[3]
    wv36f[1:cvcojmean[j]] <- as.numeric(fcst1)
    if (cvcojmean[j]<out_evhor[5]) {
      wv36f[cvcojmean[j]:out_evhor[5]]=as.numeric(fcst2)
    }
    wv33b <- wm02[j,1:(out_evhor[6] - out_evhor[7])]
    (wv36f+wv33b)
  }
}

fx_fcstgeneric <- function(fcst_mc,out_evhor,wm13){
  fcst_co <- foreach (j = 1:nrow(wm13),.combine=c("rbind")) %dopar% {
    wv36f <- rowMeans(fcst_mc[[j]][[1]])
    wv36g <- wm13[j,]
    (wv36f+wv36g)
  }
}

fx_crps_wm <- function(crps_mc,cvcojmean,out_evhor,wm02){
  crps_co <- foreach (j = 1:nrow(wm02),.combine=c("rbind")) %dopar% {
    wv36c <- matrix(ncol=out_evhor[5])
    wv36c[1:cvcojmean[j]] <- crps_mc[[j]][1,1:cvcojmean[j]]
    if (cvcojmean[j]<out_evhor[5]) {
      wv36c[cvcojmean[j]:out_evhor[5]]=crps_mc[[j]][2,cvcojmean[j]:out_evhor[5]]
    }
    wv36c
  }
}

fx_rndgrp <- function(wm01,frontierstp){
  rndgrp <- foreach (i = 1:frontierstp, .combine=c("rbind")) %:%
    foreach (j = 1:nrow(wm01), .combine=c("rbind")) %dopar%{
      rndgrp_pll <- (runif(nrow(wm01))<=(i/(frontierstp+0.3)))+0
      while (sum(rndgrp_pll)==0){rndgrp_pll <- (runif(nrow(wm01))<=(i/10))+0}
      rndgrp_pll
    }
  result <- list(rndgrp %*% wm01,rndgrp)
  return(result)
}

fx_sd_mymat <- function (mymat){
  sdev <- foreach (i = 1:nrow(mymat),.combine=c("cbind")) %do%{
    sd(mymat[i,])
  }
  return(as.numeric(sdev))
}

fx_optgrp_crps <- function (wv42){
  if (sum(wv42*wv45) > opt_min_cusd & sum(wv42*wv45) <= opt_max_cusd){
    fv01   <- rbind(wv42 %*% wm01_01 / sum(wv42),0)
    fl02   <- fx_int_fcst_kdcv(fv01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,F)
    result <- as.numeric(fl02[[2]][1,1])
  } else {result <- 10}
  return (result)
}

fx_optgrp_sdev <- function (wv42){
  if (sum(wv42*wv45) > opt_min_cusd & sum(wv42*wv45) <= opt_max_cusd){
    sd_evhor <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,crossvalsize)
    fv01     <- as.numeric(wv42 %*% wm01_01[,sd_evhor[2]:sd_evhor[1]] / sum(wv42))
    fv02     <- decompose(msts(fv01,seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
    fv03     <- fv01 - fv02$seasonal
    result   <- sd(fv03)
  } else {result <- 10}
  return (result)
}

#===========================================
# Functions Declarations: Plots
#===========================================
fx_plt_mymat <- function(wm05,myrangey){
  plot(range(1:ncol(wm05)), myrangey, bty="n", type="n")
  grid (NA,NULL, lty = 'dotted')
  par(mar=c(5,4,4,3.5), xpd=TRUE)
  colors   = rainbow(nrow(wm05))
  linetype = c(1:3)
  for (i in 1:nrow(wm05)){
    lines(wm05[i,], type="l", lwd=1.5, lty=linetype[1],col=colors[i])
  }
  # wm05mean <- colMeans(wm05)
  # lines(wm05mean, type="l", lwd=3, lty=linetype[2])
}

fx_plt_rnd_vs_opt <- function(bighlp,myrangex,myrangey,xunit) {
  plot(myrangex,myrangey, bty="n", type="n", xlab=xunit,
       ylab="Mean Demand",main=paste("optimum vs random groups for h =",bighlp[[1]]))
  grid (NA,NULL, lty = 'dotted')
  mycolors=c("darkgreen","green","darkblue","blue")
  points(bighlp[[2]][,1],bighlp[[2]][,2],col="gray80",pch=20)
  for (i in 3:length(bighlp)){
    points(bighlp[[i]][,1],bighlp[[i]][,2],col=mycolors[i-2],pch=20)
  }
  legend('topright', inset=c(0,0), legend = c("random","opt_sdev_outsample","opt_sdev_insample","opt_crps_outsample","opt_crps_insample"),
         lty=1, col=c("gray80",mycolors), bty='n', cex=.75, title="Grouping")
}

#===========================================
# Functions Declarations: Integrations
#===========================================
fx_int_fcst_kdcv <- function(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,fcst_run){
  def_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,crossvalsize)
  wm01       <- wm01_01[,def_evhor[2]:def_evhor[3]]                # work matrix
  # ------ Cross Validation ----------------
  cross_seq     <- seq(def_evhor[8],def_evhor[4]-max(ahead_t),round((def_evhor[4]-max(ahead_t)-def_evhor[7])/crossvalstps))
  crossval_runs <- foreach (k = cross_seq, .combine=c('rbind')) %do%{
    co_evhor    <- fx_evhor(wm01,k,0,ahead_t,s02,is_wins_weeks,0)
    wm01cv      <- wm01[,co_evhor[2]:co_evhor[3]]                  # work matrix
    wm02cv      <- fx_seas(wm01cv,s01,s02,sum_of_h,co_evhor)       # in-sample seasonality pattern
    wm03cv      <- wm01cv[,(co_evhor[4]+1):co_evhor[6]]            # out-sample original load data
    wm04cv      <- fx_unseas(wm01cv,wm02cv,s02,co_evhor)           # in-out sample unseasonalised
    fcst_mccv   <- fx_fcst_kds(wm04cv,win_size,co_evhor,sampling)
    crps_mccv   <- fx_crps_mc(wm04cv,fcst_mccv,co_evhor,sampling)
    cvcoj       <- fx_crossover(fcst_mccv,crps_mccv,wm02cv,co_evhor)
    wm05cv      <- fx_crps_wm(crps_mccv,cvcoj,co_evhor,wm02cv)
    c(as.numeric(cvcoj),rowMeans(wm05cv))
  }
  cvcojmean <- foreach (a = 1:nrow(wm01),.combine=c("rbind")) %dopar%{
    round(mean(crossval_runs[,a]))
  }
  cvcrpsmean <- foreach (a = (nrow(wm01)+1):(2*nrow(wm01)),.combine=c("rbind")) %dopar%{
    mean(crossval_runs[,a])
  }
  if (fcst_run == F) {
    return(list(cvcojmean,cvcrpsmean))
  }
  # ------ Forecasting & Verification ------
  if (fcst_run == T) {
    out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,0)
    wm02       <- fx_seas(wm01,s01,s02,sum_of_h,out_evhor)           # in-sample seasonality pattern
    wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]               # out-sample original load data
    wm04       <- fx_unseas(wm01,wm02,s02,out_evhor)                 # in-out sample unseasonalised
    fcst_mc    <- fx_fcst_kds(wm04,win_size,out_evhor,sampling)
    crps_mc    <- fx_crps_mc(wm04,fcst_mc,out_evhor,sampling)
    wm03fcst   <- fx_fcst_wm(fcst_mc,cvcojmean,out_evhor,wm02)
    wm05       <- fx_crps_wm(crps_mc,cvcojmean,out_evhor,wm02)
    return(list(wm03fcst,wm05,wm04[,1:out_evhor[7]]))
  }
}

fx_int_fcstgeneric_armagarch <- function(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh){
  out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,0)
  wm01       <- wm01_01[,out_evhor[2]:out_evhor[3]]                         # work matrix
  wl02       <- fx_seas2(wm01,s01,s02,sum_of_h,out_evhor)                   # in-sample seasonality pattern (s,r,t)
  wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]                        # out-sample original load data
  wm13       <- fx_unseas2(wm01,wl02,s02,out_evhor)                         # out-sample estimated trend + seas
  wm14       <- wl02[[2]]                                                   # in-sample noise
  fcst_mc    <- fx_fcst_armagarch(wm14,armalags,win_size,ahead_t,out_evhor,sampling,cross_overh) # returns list with next ahead_t fcst and sd
  wm03fcst   <- fx_fcstgeneric(fcst_mc,out_evhor,wm13)
  wm05       <- fx_crpsgeneric(wm03,wm13,wm14,fcst_mc,out_evhor,sampling)
  return(list(wm03fcst,wm05,wm14))
}

#===========================================
# BIG [h] LOOP Start
#===========================================
bighlpcrps = list()
bighlpsdev = list()

for (h in hrz_lim){
  cat("\n\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h,"\n")
  #===========================================
  # Individual customers forecast
  #===========================================
  cat("[Ind] ")
  wm01_01    <- wm01_00[min(cus_list):max(cus_list),]
  wl06kd     <- fx_int_fcst_kdcv(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T)
  wl06ag     <- fx_int_fcstgeneric_armagarch(wm01_01,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,armalags,cross_overh)
  wv45       <- rowMeans(wl06[[1]])
  sd01       <- as.numeric(fx_sd_mymat(wl06[[3]]))
  
  #===========================================
  # Random groups forecast
  #===========================================
  cat("[Rnd] ")
  wm01_02l   <- fx_rndgrp(wm01_01,frontierstp)
  wm01_02    <- wm01_02l[[1]] / rowSums(wm01_02l[[2]])
  wl06rnd    <- fx_int_fcst_kdcv(wm01_02,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T)
  wv45rnd    <- as.numeric(rowMeans(wl06rnd[[1]]) * rowSums(wm01_02l[[2]]))
  sd01rnd    <- as.numeric(fx_sd_mymat(wl06rnd[[3]]))
  cr01rnd    <- rowMeans(wl06rnd[[2]])
  
  #===========================================
  # Optimised sdev groups forecast
  #===========================================
  cat("[OptSDEV] ")
  wv46         <- seq(0,frontierstp)^2/frontierstp^2 * sum(wv45)
  optgrp_sdev  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_sdev, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             optgrp$par
                           }
  opt_min_cusd<- 0
  opt_max_cusd<- max(wv46)
  wm01_03l    <- list(optgrp_sdev %*% wm01_01, optgrp_sdev)
  wm01_03     <- wm01_03l[[1]] / rowSums(wm01_03l[[2]])
  wl06optsdev <- fx_int_fcst_kdcv(wm01_03,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T)
  wv45optsdev <- as.numeric(rowMeans(wl06optsdev[[1]]) * rowSums(wm01_03l[[2]]))
  sd01optsdev <- as.numeric(fx_sd_mymat(wl06optsdev[[1]]))
  cr01optsdev <- rowMeans(wl06optsdev[[2]])
  cr02optsdev <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud","foreach"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_crps(optgrp_sdev[i,])
                          }
  cr02optsdev <- as.numeric(cr02optsdev)
  sd02optsdev <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_sdev(optgrp_sdev[i,])
                          }
  sd02optsdev <- as.numeric(sd02optsdev)
  
  # for sdev calc: sd01 is the outsample result, sd02 is the crossval result
  # for crps calc: cr01 is the outsample result, cr02 is the crossval result
  
  #===========================================
  # Optimised crps groups forecast
  #===========================================
  cat("[OptCrossval]\n")
  optgrp_crps  <- foreach (i = 1:frontierstp,
                           .packages=c("forecast","rgenoud","foreach"),
                           .combine=c("rbind")) %dopar% {
                             opt_min_cusd  = wv46[i]
                             opt_max_cusd  = wv46[i+1]
                             optgrp   <- genoud(fx_optgrp_crps, nvars=nrow(wm01_01), max.generations=300, wait.generations=20,
                                                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                                                data.type.int=TRUE,  int.seed=1,
                                                print.level=1)
                             optgrp$par
                           }
  opt_min_cusd<- 0
  opt_max_cusd<- max(wv46)
  wm01_04l    <- list(optgrp_crps %*% wm01_01, optgrp_crps)
  wm01_04     <- wm01_04l[[1]] / rowSums(wm01_04l[[2]])
  wl06optcrps <- fx_int_fcst_kdcv(wm01_04,h,in_sample_fr,s01,s02,sum_of_h,win_size,is_wins_weeks,crossvalsize,T)
  wv45optcrps <- as.numeric(rowMeans(wl06optcrps[[1]]) * rowSums(wm01_04l[[2]]))
  sd01optcrps <- as.numeric(fx_sd_mymat(wl06optcrps[[1]]))
  cr01optcrps <- rowMeans(wl06optcrps[[2]])
  cr02optcrps <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud","foreach"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_crps(optgrp_crps[i,])
                          }
  cr02optcrps <- as.numeric(cr02optcrps)
  sd02optcrps <- foreach (i = 1:frontierstp,
                          .packages=c("forecast","rgenoud"),
                          .combine=c("cbind")) %dopar% {
                            fx_optgrp_sdev(optgrp_crps[i,])
                          }
  sd02optcrps <- as.numeric(sd02optcrps)
  
  # for sdev calc: sd01 is the outsample result, sd02 is the crossval result
  # for crps calc: cr01 is the outsample result, cr02 is the crossval result
  
  bighlpcrps[[match(h,hrz_lim)]] = list(h,cbind(cr01rnd,wv45rnd),cbind(cr01optsdev,wv45optsdev),cbind(cr02optsdev,wv45optsdev),cbind(cr01optcrps,wv45optcrps),cbind(cr02optcrps,wv45optcrps))
  bighlpsdev[[match(h,hrz_lim)]] = list(h,cbind(sd01rnd,wv45rnd),cbind(sd01optsdev,wv45optsdev),cbind(sd02optsdev,wv45optsdev),cbind(sd01optcrps,wv45optcrps),cbind(sd02optcrps,wv45optcrps))
  print(proc.time() - ptm)
}
#===========================================
# BIG [h] LOOP Plots
#===========================================
for (i in 1:length(hrz_lim)){
  fx_plt_rnd_vs_opt(bighlpcrps[[i]],c(0.02,0.07),c(0,3),"CRPS")
}

#===========================================
# Outputs
#===========================================
# saveRDS(crpsh_CLU_dyn,  file="0200_analysis.rds") # crps winsize vs time_adead, mean&desv,                per grouping 
# saveRDS(wm01_4D[[1]],   file="0200_origdat1.rds") # data customer vs time_series,                         per grouping
# saveRDS(wm02_4D[[1]],   file="0200_seasona1.rds") # seas customer vs seasonal_ts,            per horizon, per grouping
# saveRDS(fcst_4D[[1]],   file="0200_forecas1.rds") # fcst winsize vs time_ahead per customer, per horizon, per grouping
# saveRDS(parbundl0200,   file="0200_parbundl.rds")
# saveRDS(kdscrps,        file="0200_function.rds")
