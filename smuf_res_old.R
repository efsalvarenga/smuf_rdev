#===========================================
# Smart Metering Uncertainty Forecasting
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# smuf_res read and visualise results
#===========================================

setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")


#smuf_main run results

smuf_run_read  <- readRDS("smuf_run_0623_Super01.rds")
res.hlp.opgr   <- smuf_run_read[[1]]
res.hlp.crps   <- smuf_run_read[[2]]

myleg = c("Random","Sdev KDS","Sdev ARMA-GARCH","Crossval CRPS KDS","Crossval CRPS ARMA-GARCH")
for (i in 1:length(res.hlp.crps)) {
  fx_plt_rnd_vs_opt(res.hlp.crps[[i]][[2]],c(0,0.05),c(0,32),myleg,"CRPS")
}

