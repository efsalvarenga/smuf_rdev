#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM0320 optimised group forecasting
#-------------------------------------------
# Notes
# 19/04 Import optgrps
#===========================================

#===========================================
# Libraries, Inputs
#===========================================
source('./KSSM0211_read-basefcst_v01.R')
optgrp_pll   <- readRDS("0300_optgrp.rds")[[1]]
optgrp_plt   <- readRDS("0300_optgrp.rds")[[2]]
opt_ahead_t  <- readRDS("0300_optpar.rds")[1]
opt_win_sel  <- readRDS("0300_optpar.rds")[2]
opt_hrz_sel  <- readRDS("0300_optpar.rds")[3]
frontierstp  <- readRDS("0300_optpar.rds")[4]

#===========================================
# Groups Demand Calculation
#===========================================
ptm <- proc.time() # Start the clock!
cl  <- makeCluster(detectCores())
registerDoParallel(cl)

wm11 <- foreach (i = 1:frontierstp) %dopar%{
  matrix(rep((optgrp_pll[i,] %*% wm01_4D1) / sum(optgrp_pll[i,]),each=length(cus_list)), nrow=length(cus_list))
}
crpsh_CLU_opt <- kdscrps(wm11,hrz_lim[opt_hrz_sel],win_size[opt_win_sel],ahead_t[opt_ahead_t],sampling)

print(proc.time() - ptm)        # Stop the clock

#===========================================
# Outputs
#===========================================
saveRDS(crpsh_CLU_opt, file="0320_analysis.rds")
