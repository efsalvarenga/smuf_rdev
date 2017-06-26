#===========================================
# Smart Metering Uncertainty Forecasting
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#-------------------------------------------
# smuf_aux create plots for paper
#===========================================
library(ggplot2)

setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")

wm01_00       <- readRDS("smuf_import-complete.rds")
importpar     <- readRDS("smuf_import-parameter.rds")
s01           <- importpar[1]
s02           <- importpar[2]
s03           <- importpar[3]
sum_of_h      <- importpar[4]
data_size     <- importpar[5]

cus_no  <- 1000
ds_ini  <- 1001
ds_len  <- s02*4

plot1   <- as.data.frame(t(wm01_00[1:cus_no,ds_ini:(ds_ini+ds_len-1)]))
plot1   <- cbind(seq(1,ds_len),plot1)
colnames(plot1) <- c("Time",paste("Cus",1:cus_no,sep=""))
ggplot(plot1, aes(Time,Cus7)) +
  geom_line(color="firebrick") + 
  # labs(title="Example of demand pattern for one customer") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Hours (4 weeks)", y = "Energy demand (KWh)")

cus_nos <- c(3,4,8,9)
sl_win  <- seq(1993:2016)
plot2   <- as.numeric(wm14[cus_nos,sl_win])
plot2   <- as.data.frame(plot2)
plot2   <- cbind(plot2,rep(paste("Cus",cus_nos,sep=""),24))
colnames(plot2) <- c("Demand","Customer")
ggplot(plot2, aes(Demand, color=Customer, fill=Customer)) +
  geom_density(alpha = 0.1)




saveRDS(list(plot1,plot2),file="smuf_aux_plots")
