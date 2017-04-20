#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM0330 opt group vs random
#-------------------------------------------
# Notes
# 11/04 Modular implemented
# 18/04 Parallel processing and data export
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
library(doParallel)
library(rgenoud)

ptm <- proc.time() # Start the clock!
cl  <- makeCluster(detectCores(),outfile="")
registerDoParallel(cl)

source('./KSSM0211_read-basefcst_v01.R')

#===========================================
# Desvcalc functions
#===========================================
cus_grouping_B <- function(wv42){
  fv01   <- (wv42 %*% (wm01_4D1[,(in_sample_ini):(event_horizon)])) / sum(wv42)
  fv02   <- decompose(msts(fv01[1,],seasonal.periods=c(s01/sum_of_h,s02/sum_of_h)))
  fv03   <- fv02$random
  randev <- sd(fv03,na.rm=TRUE)
}

#===========================================
# Random Grouping
#===========================================
opt_ahead_t   = 1
opt_win_sel   = 1
opt_hrz_sel   = 2
frontierstp   = 20
rndsim_size   = 100

rndgrp = matrix(nrow=rndsim_size,ncol=3)
rndgrp_plt <- foreach (i = 1:10, .packages=c("forecast"), .combine=c("rbind")) %dopar%{
  for (j in 1:rndsim_size){
    rndgrp_pll = (runif(length(cus_list))<=(i/10))+0
    rndgrp[j,1] = sum(rndgrp_pll*wv45)
    rndgrp[j,2] = cus_grouping_B(rndgrp_pll)
    rndgrp[j,3] = sum(rndgrp_pll)
  }
  rndgrp
}
rndgrp_pltres = matrix(nrow=frontierstp,ncol=3)
wv46          = seq(0,sum(wv45),sum(wv45)/frontierstp)

for (i in 1:frontierstp){
  print(i)
  rndgrp_pltres[i,] = colMeans(rndgrp_plt[(rndgrp_plt[,1]>wv46[i]&rndgrp_plt[,1]<wv46[i+1]),])
}

plot(range(1:5/100),range(0:round(sum(wv45))), bty="n", type="n", xlab=paste("stdev of noise after decompose"),
     ylab="Mean Demand",main=paste("optimum vs random groups"))
grid (NA,NULL, lty = 'dotted')
points(optgrp_plt[2:20,2],optgrp_plt[2:20,1],col="green",pch=19)
points(rndgrp_pltres[2:20,2],rndgrp_pltres[2:20,1],col="red",pch=19)
legend('topright', inset=c(-0.10,0), legend = c("random","optimum"),
       lty=1, col=c("red","green"), bty='n', cex=.75, title="Groups")
