#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM9100 plots
#-------------------------------------------
# Notes
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
source('./KSSM0211_read-basefcst_v01.R')
source('./KSSM0321_read-optgfcst_v01.R')

#===========================================
# Plot Parameters
#===========================================
sim_soh       = sum_of_h
sim_ahead_t_n = max(ahead_t)
sim_rangefcst = sim_ahead_t_n/sim_soh
sim_cus_clu   = c("ind",cus_clu)
sim_winsiz_n  = 3
range_ch1     = range(0,0.12)

#===========================================
# Organising imported data
#===========================================
sim_rangefcst = sim_ahead_t_n/sim_soh
main_lists    = length(crpsh_CLU_dyn)
wm21_automean = matrix(nrow=main_lists,ncol=sim_ahead_t_n)
wm21_autodesv = matrix(nrow=main_lists,ncol=sim_ahead_t_n)
for (i in 1:main_lists) {
  wm21_automean[i,] = crpsh_CLU_dyn[[i]][[1]][sim_winsiz_n,]
  wm21_autodesv[i,] = crpsh_CLU_dyn[[i]][[2]][sim_winsiz_n,]
}
optg_lists    = length(crpsh_CLU_opt)
wm22_optgmean = matrix(nrow=optg_lists,ncol=5)
wm22_optgdesv = matrix(nrow=optg_lists,ncol=5)
for (i in 1:optg_lists) {
  wm22_optgmean[i,1] = crpsh_CLU_opt[[i]][[1]][1]
  wm22_optgmean[i,2] = opt_ahead_t
  wm22_optgmean[i,3] = optgrp_plt[i,1]
  wm22_optgmean[i,4] = optgrp_plt[i,2]
  wm22_optgmean[i,5] = optgrp_plt[i,3]
}

#===========================================
# Plot of wm21_automean
#===========================================
plot(range(1:(sim_rangefcst)), range_ch1, bty="n", type="n", xlab=paste("Time-blocks ahead forecast"),
     ylab="CRPS mean",main=paste("CRPS for Groups, window size 4h to 24h"))
grid (NA,NULL, lty = 'dotted')
par(mar=c(5,4,4,3.5), xpd=TRUE)
colors <- rainbow(length(sim_cus_clu))
linetype <- c(1:3)
for (i in 1:(length(sim_cus_clu))){
  lines(wm21_automean[i,(1:sim_rangefcst)], type="l", lwd=1.5, lty=linetype[1],col=colors[i])
}
# points(wm22_optgmean[,2],wm22_optgmean[,1],col="red",pch=19)
# legend('topright', inset=c(-0.10,0), legend = c(sim_cus_clu),
#        lty=1, col=rainbow(length(sim_cus_clu)), bty='n', cex=.75, title="Groups")
# 
# plot(wm22_optgmean[,5],wm22_optgmean[,1],main=paste("CRPS vs #customer in group"))
# plot(wm22_optgmean[,5],wm22_optgmean[,4],main=paste("sdev vs #customer in group"))
# plot(wm22_optgmean[,3],wm22_optgmean[,1],main=paste("CRPS vs aggr demand /group"))
# plot(wm22_optgmean[,3],wm22_optgmean[,4],main=paste("sdev vs aggr demand /group"))
# plot(wm22_optgmean[,4],wm22_optgmean[,1],main=paste("sdev vs CRPS"))
