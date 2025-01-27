#===========================================
# Smart Metering Uncertainty Forecasting
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# smuf_lags checking lags for ARMA-GARCH
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
library(forecast)
library(doParallel)
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
cus_list      <- seq(1,1000)
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(0,9)*2069
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
maxlag        <- 10                      # Max lags analysed for ARIMA fit (ARMA-GARCH model)

#===========================================
# Functions Declarations
#===========================================
fx_lags_armagarch <- function (wm14,maxlag,out_evhor){
  lags_armagarch  <- foreach (j = 1:nrow(wm14), .packages=c("rugarch"), .combine=c("rbind")) %dopar% {
    runvec        <- wm14[j,]
    # Defining ARMA lags
    final.bic <- matrix(nrow=0,ncol=4)
    for (p in 0:maxlag) for (q in 0:maxlag) {
      if ( p == 0 && q == 0) {
        next
      }
      arimaFit = tryCatch(arima(runvec, order=c(p, 0, q)),
                          error=function( err ) FALSE,
                          warning=function( err ) FALSE )
      if( !is.logical( arimaFit ) ) {
        final.bic <- rbind(final.bic,c(p,q,AIC(arimaFit,k=log(out_evhor[7])),AIC(arimaFit)))
      } else {
        next
      }
    }
    final.ord <- final.bic[sort.list(final.bic[,3]), ]
    if (nrow(final.ord)==0) {
      armapar <- c(1,1)
    } else {
      armapar <- c(final.ord[1,1],final.ord[1,2])
    }
    armapar
  }
  return(lags_armagarch)
}

#===========================================
# BIG [h] LOOP Start
#===========================================
wm01_01    <- wm01_00[min(cus_list):length(cus_list),]
hlags      <- matrix(nrow=0,ncol=2)
for (h in hrz_lim) {
  cl  <- makeCluster(detectCores())   # reset parallel workers
  registerDoParallel(cl)
  cat("\nStep",match(h,hrz_lim), "of",length(hrz_lim),"| Running BIG [h] LOOP with h =",h)
  out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,0)
  wm01       <- wm01_01[,out_evhor[2]:out_evhor[3]]                       # work matrix
  # wm01       <- wm01[-(which(!is.na(match(rowMeans(wm01), 0)))),]       # remove customers with no data (not in use)  
  wl02       <- fx_seas2(wm01,s01,s02,sum_of_h,out_evhor)                   # in-sample seasonality pattern (s,r,t)
  wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]                        # out-sample original load data
  wm13       <- fx_unseas2(wm01,wl02,s02,out_evhor)                         # out-sample estimated trend + seas
  wm14       <- wl02[[2]]                                                   # in-sample noise
  hlag       <- fx_lags_armagarch(wm14,maxlag,out_evhor)
  hlags      <- rbind(hlags,hlag)
  # cat("\nlags",hlag[1],hlag[2],"customers",hlag[3],"\n")
  saveRDS(hlags,  file="smuf_lags-arma.rds")
  print(proc.time() - ptm)
}

#===========================================
# Outputs
#===========================================
saveRDS(hlags,  file="smuf_lags-arma.rds")
