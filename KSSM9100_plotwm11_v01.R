#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#===========================================
# Xmod9100 plot of wm11 result matrix
#===========================================
# Notes
# Modular implemented in 11/Apr/17
#
# Output
# Plots of wm result matrix
#===========================================

# source('~/_mysync/load_forecasting/KSSM_Xmod0210_summary-basefcst_v01.R')

#===========================================
# Parameters & Start data
#===========================================
sim_soh       = sum_of_h
sim_ahead_t_n = max(ahead_t)
sim_rangefcst = sim_ahead_t_n/sim_soh
sim_cus_clu   = c("ind",cus_clu)
sim_winsiz_n  = 3
range_ch1     = range(0,0.12)
range_ch2     = range(0.08,0.14)

#===========================================
# Plot of wm11_automean
#===========================================
plot(range(1:(sim_rangefcst)), range_ch1, bty="n", type="n", xlab=paste("Time-blocks ahead forecast"),
     ylab="CRPS mean",main=paste("CRPS for Groups, window size 4h to 24h"))
# axis(1, at=(1:(length(win_size))), labels=win_size)
grid (NA,NULL, lty = 'dotted')
par(mar=c(5,4,4,3.5), xpd=TRUE)
colors <- rainbow(length(sim_cus_clu))
linetype <- c(1:3)
for (i in 1:(length(sim_cus_clu))){
  lines(wm11_automean[i,(1:sim_rangefcst)], type="l", lwd=1.5, lty=linetype[1],col=colors[i])
}
legend('topright', inset=c(-0.10,0), legend = c(sim_cus_clu),
       lty=1, col=rainbow(length(sim_cus_clu)), bty='n', cex=.75, title="Groups")