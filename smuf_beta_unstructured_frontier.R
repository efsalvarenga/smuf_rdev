a <- foreach (i=1:100, .combine='rbind', .packages=c("doParallel")) %dopar% {
  binvec <- round(runif(20))
  mycrps <- fx_optgrp_crps(binvec,F)
  grpdem <- sum(binvec*wv45)
  c(binvec,mycrps,grpdem)
}



max.gen       <- 800                     # For genetic opt
waitgen       <- 500                      # For genetic opt

optuf <- genoud(fx_optgrp_sdev, nvars=nrow(wm01_01), max.generations=max.gen, wait.generations=waitgen,
                starting.values=c(rep(1,nrow(wm01_01))), Domains = cbind(c(rep(0,nrow(wm01_01))),c(rep(1,nrow(wm01_01)))),
                data.type.int=TRUE,  int.seed=1,
                print.level=1)

a <- fx_optgrp_crps(optuf$par,F)
b <- sum(optuf$par*wv45)
