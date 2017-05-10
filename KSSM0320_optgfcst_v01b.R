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
optgrp_pllA  <- readRDS("0300_optgrp_sdev.rds")[[1]]
optgrp_pltA  <- readRDS("0300_optgrp_sdev.rds")[[2]]
optgrp_pllC  <- readRDS("0300_optgrp_crps.rds")[[1]]
optgrp_pltC  <- readRDS("0300_optgrp_crps.rds")[[2]]
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

wm11A <- foreach (i = 1:frontierstp) %dopar%{
  matrix(rep((optgrp_pllA[i,] %*% wm01_4D1) / sum(optgrp_pllA[i,]),each=length(cus_list)), nrow=length(cus_list))
}
wm11C <- foreach (i = 1:frontierstp) %dopar%{
  matrix(rep((optgrp_pllC[i,] %*% wm01_4D1) / sum(optgrp_pllC[i,]),each=length(cus_list)), nrow=length(cus_list))
}

crpsh_CLU_optA <- kdscrps(wm11A,hrz_lim[opt_hrz_sel],win_size[opt_win_sel],ahead_t[opt_ahead_t],in_sample_fr, crossvalsize, sampling)
crpsh_CLU_optC <- kdscrps(wm11C,hrz_lim[opt_hrz_sel],win_size[opt_win_sel],ahead_t[opt_ahead_t],in_sample_fr, crossvalsize, sampling)

crpsh_CLU_optA[[1]][[1]]
crpsh_CLU_optC[[1]][[1]]

crpsh_CLU_optA[[2]][[1]]
crpsh_CLU_optC[[2]][[1]]

crpsh_CLU_optA[[3]][[1]]
crpsh_CLU_optC[[3]][[1]]

crpsh_CLU_optA[[4]][[1]]
crpsh_CLU_optC[[4]][[1]]

print(proc.time() - ptm)        # Stop the clock

#===========================================
# Outputs
#===========================================
saveRDS(crpsh_CLU_opt, file="0320_analysis.rds")
