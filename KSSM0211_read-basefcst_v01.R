# teste comments

#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# KSSM0211 reading base forecasting
#-------------------------------------------
# Notes
#===========================================
library(forecast)
library(verification)
library(doParallel)

crpsh_CLU_dyn  <- readRDS("0200_analysis.rds")
wm01_4D1       <- readRDS("0200_origdat1.rds")
wm02_4D1       <- readRDS("0200_seasona1.rds")
fcst_4D1       <- readRDS("0200_forecas1.rds")
parbundl0200   <- readRDS("0200_parbundl.rds")
kdscrps        <- readRDS("0200_function.rds")

s01            <- parbundl0200[[1]][1]
s02            <- parbundl0200[[1]][2]
s03            <- parbundl0200[[1]][3]
sum_of_h       <- parbundl0200[[1]][4]
data_size      <- parbundl0200[[1]][5]
in_sample_fr   <- parbundl0200[[1]][6]
cus_list       <- parbundl0200[[2]]
cus_clu        <- parbundl0200[[3]]
win_size       <- parbundl0200[[4]]
ahead_t        <- parbundl0200[[5]]
hrz_lim        <- parbundl0200[[6]]
sampling       <- parbundl0200[[7]]

