ptm <- proc.time() # Start the clock!

cus_list      <- seq(1,20)
wm01_01    <- wm01_00[min(cus_list):length(cus_list),]
h = 5000

wm01       <- wm01_01[,def_evhor[2]:def_evhor[3]]                # work matrix
out_evhor  <- fx_evhor(wm01_01,h,in_sample_fr,s02,seas_bloc_ws,0)
wm02       <- fx_seas(wm01,s01,s02,sum_of_h,out_evhor)           # in-sample seasonality pattern
wm03       <- wm01[,(out_evhor[4]+1):out_evhor[6]]               # out-sample original load data
wm04       <- fx_unseas(wm01,wm02,s02,out_evhor)                 # in-out sample unseasonalised
fcst_mc    <- fx_fcst_kds(wm04,win_size,out_evhor,sampling)

testecomb <- foreach (j = 1:nrow(wm04), .combine=c("rbind"), .packages=c("verification","doParallel","rugarch")) %dopar% {
  testec = j
  testvec = wm04[testec,]
  testkdfcst = fcst_mc[[testec]][,3:4]
  
  teste_is = testvec[1:out_evhor[7]]
  teste_os = testvec[(out_evhor[7]+1):out_evhor[6]]
  
  final.aic <- Inf
  final.order <- c(0,0,0)
  for (p in 0:5) for (q in 0:5) {
    if ( p == 0 && q == 0) {
      next
    }
    arimaFit = tryCatch(arima(teste_is, order=c(p, 0, q)),
                        error=function( err ) FALSE,
                        warning=function( err ) FALSE )
    if( !is.logical( arimaFit ) ) {
      current.aic <- AIC(arimaFit)
      if (current.aic < final.aic) {
        final.aic <- current.aic
        final.order <- c(p, 0, q)
        final.arima <- arima(teste_is, order=final.order)
      }
    } else {
      next
    }
  }
  
  # setting GARCH spec (1,1)
  spec = ugarchspec(variance.model=list(garchOrder=c(1,1)),
                    mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
                    distribution.model="sged")
  
  # fitting ARMA-GARCH parameters
  fit = tryCatch(ugarchfit(spec, teste_is, solver = 'hybrid'),
                 error=function(e) e, warning=function(w) w)
  
  sim1 = ugarchsim(fit, n.sim = win_size[2], m.sim = sampling)                # density simulation forecast
  
  testagfcst = list()
  for (i in 1: win_size[2]) {
    testagfcst[[i]] = data.frame(sim1@simulation$seriesSim[i,],sim1@simulation$sigmaSim[i,])
  }
  
  testkdcrps <- foreach (i = 1:out_evhor[5], .combine=c("cbind"), .packages=c("verification")) %dopar% {
    crps(rep(teste_os[i],sampling),testkdfcst)$CRPS
  }
  testagcrps <- foreach (i = 1:out_evhor[5], .combine=c("cbind"), .packages=c("verification")) %dopar% {
    crps(rep(teste_os[i],sampling),testagfcst[[1]])$CRPS
  }
  
  c(mean(testkdcrps),mean(testagcrps))
}

print(proc.time() - ptm)

# based on: https://www.quantstart.com/articles/ARIMA-GARCH-Trading-Strategy-on-the-SP500-Stock-Market-Index-Using-R