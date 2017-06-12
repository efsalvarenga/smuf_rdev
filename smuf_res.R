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

smuf_run_read  <- readRDS("smuf_run_0611_1400.rds")
res.bighlpopgr <- smuf_run_read[[1]]
res.bighlpcrps <- smuf_run_read[[2]]
