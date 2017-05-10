kdscrps <- function(wm01_4D, hrz_lim, win_size, ahead_t, in_sample_fr, crossvalsize, sampling){
  c_group     = seq(1,length(wm01_4D))
  crpsh_CLU   = list()
  wm02_4D   <<- list()
  fcst_4D   <<- list()
  for (c in c_group){
    cat("\n\n\n----------------------------------------\nGROUP NUMBER ",c," OF ",max(c_group),"\n----------------------------------------")
    wm01    =  wm01_4D[[c]]
    fcsth   =  list()
    crpsh   =  list()
    wm02_3D =  list()
    for (h in hrz_lim) {
      cat("\nEvent Horizon +",h,"for",nrow(wm01),"customers: ")
      # WorkMatrix 02, 03 & 04 ===================
      # wm02 in-sample seasonality pattern
      # wm03 out-sample load data
      # wm04 in- and out-sample de-seasonalised
      event_horizon = ncol(wm01)*in_sample_fr+1 + h -s02
      in_sample_ini = event_horizon - min(seas_bloc_ws,(event_horizon %/% s02)) * s02 + 1
      outsample_end = (ncol(wm01)*(1-in_sample_fr)) %/% s02 * s02 + event_horizon - (h%/%s02*s02)
      in_sample_siz = event_horizon - in_sample_ini + 1
      outsample_siz = outsample_end - event_horizon
      inosample_siz = in_sample_siz + outsample_siz
      event_hrz_mod = event_horizon - in_sample_ini + 1
      wm02          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=s02)
      wm03          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=outsample_siz)
      wm04          = matrix(nrow=nrow(wm01_4D[[c]]), ncol=(inosample_siz))
      for (j in 1:nrow(wm01_4D[[c]])){
        if (j == 1){cat("[Seasonality] ")}
        # cat(".")
        wv31i  = wm01[j,(in_sample_ini):(event_horizon)]
        wv31o  = wm01[j,(event_horizon+1):(outsample_end)]
        wv32i  = decompose(msts(wv31i,seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
        wv32is = wv32i$seasonal[1:(s02)]  #####shoud included the /sum_of_h, and 3 lines below?
        wm02[j,1:s02] = wv32is
        wm03[j,1:outsample_siz] = wv31o
        wm04[j,1:inosample_siz] = c(wv31i,wv31o) - rep(wv32is,(inosample_siz)/s02)
      }
      wm02_3D[[match(h,hrz_lim)]] = wm02
      # Forecasting and CRPS =====================
      cat("[Forecasting]")
      fcstcrpsj <- foreach (j = 1:nrow(wm01_4D[[c]]),
                            .packages=c("forecast","verification")) %dopar% {
                              # cat(".")
                              crpski = matrix(nrow=length(win_size),ncol=length(ahead_t))
                              fcstki = matrix(nrow=length(win_size),ncol=length(ahead_t))
                              for (k in win_size){
                                for (i in ahead_t){
                                  densi = density(wm04[j,(event_hrz_mod - k + 1):(event_hrz_mod)])
                                  fcsti = sample(densi$x, sampling, replace=TRUE, prob=densi$y)
                                  crpsi = crps(rep(wm04[j,(event_hrz_mod + i)], sampling), data.frame(fcsti, densi$bw))
                                  crpski[match(k,win_size),i] = crpsi$CRPS
                                  fcstki[match(k,win_size),i] = mean(fcsti)
                                }
                              }
                              list(fcstki,crpski)
                            }
      fcstkij = list()
      crpskij = list()
      for (j in 1:nrow(wm01_4D[[c]])){
        fcstkij[[j]] = fcstcrpsj[[j]][[1]]
        crpskij[[j]] = fcstcrpsj[[j]][[2]]
      }
      crpsh[[match(h,hrz_lim)]]  =  crpskij
      fcsth[[match(h,hrz_lim)]]  =  fcstkij
    }
    fcst_4D[[c]] <<- fcsth
    wm02_4D[[c]] <<- wm02_3D
    # Statistical data from CRPS 4d results =========
    crpsh_All = list()
    crpsh_All_meanki = matrix(nrow=length(win_size),ncol=length(ahead_t))
    crpsh_All_desvki = matrix(nrow=length(win_size),ncol=length(ahead_t))
    cat('\nCreating statistical data, removing dimention of Horizon and Customers')
    for (k in win_size){
      for (i in ahead_t){
        crpsh_All_tempki = integer(nrow(wm01_4D[[c]])*length(hrz_lim))
        for (h in hrz_lim){
          for (j in 1:nrow(wm01_4D[[c]])){
            crpsh_All_tempki[j+(nrow(wm01_4D[[c]]))*(match(h,hrz_lim)-1)] = crpsh[[match(h,hrz_lim)]][[j]][match(k,win_size),i]
          }
        }
        crpsh_All_meanki[match(k,win_size),i] = mean(crpsh_All_tempki)
        crpsh_All_desvki[match(k,win_size),i] = sd(crpsh_All_tempki)
      }
    }
    crpsh_All = list(crpsh_All_meanki,crpsh_All_desvki)
    crpsh_CLU[[c]] = crpsh_All
  }
  return(crpsh_CLU)
}
