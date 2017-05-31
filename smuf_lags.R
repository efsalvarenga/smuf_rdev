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
cus_list      <- seq(1,12)
ahead_t       <- seq(1, (24/sum_of_h))   # Up to s02
hrz_lim       <- seq(0,1)*2069
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
seas_bloc_ws  <- 6                       # Number of weeks used for calculating seasonality pattern (6 seems best)
maxlag        <- 7                       # Max lags analysed for ARIMA fit (ARMA-GARCH model)

#===========================================
# Functions Declarations
#===========================================
fx_lags_armagarch <- function (wm04,maxlag,out_evhor){
  lags_armagarch <- foreach (j = 1:nrow(wm04), .packages=c("rugarch"), .combine=c("rbind")) %dopar% {
    runvec       <- wm04[j,1:out_evhor[7]]
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
    c(final.ord[1,1],final.ord[1,2])
  }
  return(c(max(lags_armagarch[,1]),max(lags_armagarch[,2]),nrow(lags_armagarch)))
}

#===========================================
# BIG [h] LOOP Start
#===========================================
wm01_01    <- wm01_00[min(cus_list):length(cus_list),]
hlags      <- foreach (h = hrz_lim, .packages=c("rugarch","doParallel"), .combine=c("rbind")) %dopar% {
  out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,seas_bloc_ws,0)
  wm01       <- wm01_01[,out_evhor[2]:out_evhor[3]]                       # work matrix
  wm02       <- fx_seas(wm01,s01,s02,sum_of_h,out_evhor)                  # in-sample seasonality pattern
  wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]                      # out-sample original load data
  wm04       <- fx_unseas(wm01,wm02,s02,out_evhor)                        # in-out sample unseasonalised
  fx_lags_armagarch(wm04,maxlag,out_evhor)
}
arma.lags <- c(max(hlags[,1]),max(hlags[,2]))
cus_check <- hlags[,3]

#===========================================
# Outputs
#===========================================
saveRDS(arma.lags,  file="smuf_lags-arma.rds")
print(proc.time() - ptm)
