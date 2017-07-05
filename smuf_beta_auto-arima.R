# smuf_beta_auto.arima

closeAllConnections()
setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")
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
# frontierstp   <- stepwiseHPC             # Number of demand bins (Stepwise frontier for portfolio optimisation)
# frontierexp   <- 1.2                     # Exponentiality of frontier steps
# max.gen       <- 300                     # For genetic opt
# waitgen       <- 50                      # For genetic opt
win_size      <- c(4,24)                 # Small and large win_size (select only 2)
win_selec     <- win_size[1]
# cross_overh   <- 4                       # Cross-over forced for fx_fcst_kds_quickvector
ahead_t       <- seq(1, (12/sum_of_h))    # Up to s02
h             <- 0
hrz_lim       <- seq(5,6)*113            # Rolling forecasts steps {seq(0:167)*113} is comprehensive
in_sample_fr  <- 1/6                     # Fraction for diving in- and out-sample
# crossvalsize  <- 1                       # Number of weeks in the end of in_sample used for crossvalidation
# crossvalstps  <- 16                      # Steps used for multiple crossvalidation (Only KDE)
# crossvalfocus <- c(1,2)                  # What period is focused when running crossvalidation
is_wins_weeks <- 12                      # Number of weeks used for in-sample (KDE uses win_size) & seasonality
sampling      <- 1024                    # For monte-carlo CRPS calculation
armalags      <- c(3,3)                  # Max lags for ARIMA fit in ARMA-GARCH model (use smuf_lags.R)
gof.min       <- 0                    # GoF crossover value to change ARMA-GARCH to KDS

wm01_01    <- wm01_00[min(cus_list):max(cus_list),]

out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,ahead_t,s02,is_wins_weeks,0)
wm01       <- wm01_01[,out_evhor[2]:out_evhor[3]]                         # work matrix
wl02       <- fx_seas2(wm01,s01,s02,sum_of_h,out_evhor)                   # in-sample seasonality pattern (s,r,t)
wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]                        # out-sample original load data
wm13       <- fx_unseas2(wm01,wl02,s02,out_evhor)                         # out-sample estimated trend + seas
wm14       <- wl02[[2]]                                                   # in-sample noise

# compare arma_lags time
# ptm    <- proc.time()
# efsa_armalags  <- foreach (j = 1:nrow(wm14), .packages=c("rugarch","doParallel"), .combine='rbind') %dopar% {
#   runvec       <- wm14[j,]
#   # Defining ARMA lags
#   final.bic <- matrix(nrow=0,ncol=4)
#   for (p in 0:armalags[1]) for (q in 0:armalags[2]) {
#     if ( p == 0 && q == 0) {
#       next
#     }
#     arimaFit = tryCatch(arima(runvec, order=c(p, 0, q)),
#                         error=function(err) FALSE,
#                         warning=function( err ) FALSE )
#     if( !is.logical(arimaFit) ) {
#       final.bic <- rbind(final.bic,c(p,q,AIC(arimaFit,k=log(out_evhor[7])),AIC(arimaFit)))
#     } else {
#       next
#     }
#   }
#   final.ord <- final.bic[sort.list(final.bic[,3]), ]
#   c(final.ord[1,1],0,final.ord[1,2:3])
# }
# print(proc.time() - ptm)
# 
# ptm    <- proc.time()
# auto_armalags  <- foreach (j = 1:nrow(wm14), .packages=c("rugarch","doParallel"), .combine='rbind') %dopar% {
#   runvec       <- wm14[j,]
#   autotest     <- auto.arima(runvec,d=0,D=0,max.p=armalags[1],max.q=armalags[2],ic='bic')
#   c(arimaorder(autotest),autotest$bic)
# }
# print(proc.time() - ptm)

#compare CRPS
# ptm    <- proc.time()
# fcst_mc    <- fx_fcst_armagarch(wm14,armalags,win_selec,ahead_t,out_evhor,sampling,cross_overh,gof.min) # returns list with next ahead_t fcst and sd
# wm03fcst   <- fx_fcstgeneric(fcst_mc,out_evhor,wm13)
# wm05       <- fx_crpsgeneric(wm03,wm13,wm14,fcst_mc,out_evhor,sampling)
# print(proc.time() - ptm)

ptm    <- proc.time()
fcst_mcs    <- fx_fcst_armagarchs(wm14,(armalags*3),win_selec,ahead_t,out_evhor,sampling,cross_overh,gof.min) # returns list with next ahead_t fcst and sd
wm03fcsts   <- fx_fcstgeneric(fcst_mcs,out_evhor,wm13)
wm05s       <- fx_crpsgeneric(wm03,wm13,wm14,fcst_mcs,out_evhor,sampling)
print(proc.time() - ptm)

ptm    <- proc.time()
fcst_mca    <- fx_fcst_armagarchs(wm14,(armalags),win_selec,ahead_t,out_evhor,sampling,cross_overh,gof.min) # returns list with next ahead_t fcst and sd
wm03fcsta   <- fx_fcstgeneric(fcst_mca,out_evhor,wm13)
wm05a       <- fx_crpsgeneric(wm03,wm13,wm14,fcst_mca,out_evhor,sampling)
print(proc.time() - ptm)



colMeans(wm05a) - colMeans(wm05s)

