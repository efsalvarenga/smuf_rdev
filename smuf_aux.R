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

plt1nam <- "compare_demands.pdf"
cus_no  <- 1000
ds_ini  <- 1001
ds_len  <- s02*4
plot1   <- as.data.frame(t(wm01_00[1:cus_no,ds_ini:(ds_ini+ds_len-1)]))
plot1   <- cbind(seq(1,ds_len),plot1)
colnames(plot1) <- c("Time",paste("Cus",1:cus_no,sep=""))
ggplot1a <- ggplot(plot1, aes(Time,Cus2)) +
  geom_line(color="firebrick") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Hours (4 weeks)", y = "Energy demand (KWh)") +
  theme(text=element_text(family="Times"))
ggplot1b <- ggplot(plot1, aes(Time,Cus4)) +
  geom_line(color="firebrick") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Hours (4 weeks)", y = "Energy demand (KWh)") +
  theme(text=element_text(family="Times"))
ggplot1c <- ggplot(plot1, aes(Time,Cus5)) +
  geom_line(color="firebrick") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Hours (4 weeks)", y = "Energy demand (KWh)") +
  theme(text=element_text(family="Times"))
ggplot1d <- ggplot(plot1, aes(Time,Cus6)) +
  geom_line(color="firebrick") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Hours (4 weeks)", y = "Energy demand (KWh)") +
  theme(text=element_text(family="Times"))
multiplot(ggplot1a, ggplot1d, ggplot1c, ggplot1b, cols=2)
# ggsave(paste(Sys.Date(),plt1nam,sep="_"),path="./Plots")

plt2nam <- "compare_densities.pdf"
cus_nos <- c(3,4,8,9)
sl_win  <- seq(1993:2016)
plot2   <- as.numeric(wm14[cus_nos,sl_win])
plot2   <- as.data.frame(plot2)
plot2   <- cbind(plot2,rep(paste("Cus",cus_nos,sep=""),24))
colnames(plot2) <- c("Demand","Customer")
ggplot2 <- ggplot(plot2, aes(Demand, color=Customer, fill=Customer)) +
            geom_density(alpha = 0.1) +
            theme(text=element_text(family="Times"))
ggplot2
# ggsave(paste(Sys.Date(),plt2nam,sep="_"),path="./Plots")

plt3nam <- "optgrp_01.pdf"
plot3i  <- readRDS("smuf_run_0624_defheurd.rds")
myleg   <- c("Random","SDKD","SDAG","CVKD","CVAG")
fx_plt_rnd_vs_opt(plot3i[[2]][[length(plot3i[[2]])]][[2]],c(0,0.15),c(0,15),myleg,"CRPS")
plot3   <- cbind(as.data.frame(plot3i[[2]][[1]][[2]][[2]]),myleg[1])
colnames(plot3) <- c("CRPS","uDemand","Grouping")
niceleg <- c("Optimal (st. dev.)","","Optimal (cross-val.)")
for (i in c(4,6)) {
  temp  <- cbind(as.data.frame(plot3i[[2]][[1]][[2]][[i]]),niceleg[(i-3)])
  colnames(temp) <- c("CRPS","uDemand","Grouping")
  plot3 <- rbind(plot3,temp)
}
ggplot3 <- ggplot(plot3, aes(CRPS,uDemand, color=Grouping)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray90"),
                     panel.grid.minor = element_line(colour = "gray95"), axis.line = element_line(colour = "gray60")) +
  scale_color_manual(values=c("gray80", "dodgerblue3", "firebrick")) +
  theme(text=element_text(family="Times",size=18)) +
  scale_y_continuous(name="Mean Demand (in KWh)") +
  scale_x_continuous(name="Forecast Uncertainty (in average kWh CRPS)",
                     limits=c(0, 0.1),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))
ggplot3
ggsave(paste(Sys.Date(),plt3nam,sep="_"),path="./Plots")

#===========================================
saveRDS(list(plot1,plot2,plot3),file="smuf_aux_plots.rds")