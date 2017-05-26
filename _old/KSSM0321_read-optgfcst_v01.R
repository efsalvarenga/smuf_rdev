#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM0321 reading optgrp forecasting
#-------------------------------------------
# Notes
#===========================================
optgrp_pll    <- readRDS("0300_optgrp.rds")[[1]]
optgrp_plt    <- readRDS("0300_optgrp.rds")[[2]]
opt_ahead_t   <- readRDS("0300_optpar.rds")[1]
opt_win_sel   <- readRDS("0300_optpar.rds")[2]
opt_hrz_sel   <- readRDS("0300_optpar.rds")[3]
frontierstp   <- readRDS("0300_optpar.rds")[4]
crpsh_CLU_opt <- readRDS("0320_analysis.rds")
