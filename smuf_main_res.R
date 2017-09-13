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
library(reshape)
library(dplyr)
library(magrittr)
library(data.table)

setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")

wm01_00       <- readRDS("smuf_import-complete.rds")
importpar     <- readRDS("smuf_import-parameter.rds")
implotdt      <- readRDS("smuf_aux_plots.rds")
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
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),  
        axis.title.x = element_text(color="black",size=18),
        axis.title.y = element_text(color="black",size=18))
ggplot1b <- ggplot(plot1, aes(Time,Cus4)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),  
        axis.title.x = element_text(color="black",size=18),
        axis.title.y = element_text(color="black",size=18))
ggplot1c <- ggplot(plot1, aes(Time,Cus5)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),  
        axis.title.x = element_text(color="black",size=18),
        axis.title.y = element_text(color="black",size=18))
ggplot1d <- ggplot(plot1, aes(Time,Cus6)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),  
        axis.title.x = element_text(color="black",size=18),
        axis.title.y = element_text(color="black",size=18))
multiplot(ggplot1a, ggplot1d, ggplot1c, ggplot1b, cols=2)
# ggsave(paste(Sys.Date(),plt1nam,sep="_"),path="./Plots")

plt2nam <- "compare_densities.pdf"
plot2   <- implotdt[[2]]
# cus_nos <- c(3,4,8,9)
# sl_win  <- seq(1993:2016)
# plot2   <- as.numeric(wm14[cus_nos,sl_win])
# plot2   <- as.data.frame(plot2)
# plot2   <- cbind(plot2,rep(paste("Cus",cus_nos,sep=""),24))
# colnames(plot2) <- c("Demand","Customer")
ggplot2 <- ggplot(plot2, aes(Demand, linetype=Customer)) +
            theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
            geom_density(alpha = 0.1) +
            theme(text=element_text(family="Times"),
                  axis.text.x = element_text(color="black",size=18),
                  axis.text.y = element_text(color="black",size=18),  
                  axis.title.x = element_text(color="black",size=18),
                  axis.title.y = element_text(color="black",size=18)) +
            theme(legend.position="none") +
            labs(x = "Deseasonalised Demand (kWh)", y = "Density") +
            scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0),limits = c(0,12.5))
ggplot2
# ggsave(paste(Sys.Date(),plt2nam,sep="_"),path="./Plots")

plt3nam <- "optgrp_01.pdf"
# plot3i  <- readRDS("smuf_run_0624_defheurd.rds")
plot3i  <- readRDS("smuf_run_0704_defheur_med01.rds")
myleg   <- c("Random","SDKD","SDAG","CVKD","CVAG")
fx_plt_rnd_vs_opt(plot3i[[2]][[length(plot3i[[2]])]][[2]],c(0,0.05),c(0,27),myleg,"CRPS")
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
  scale_y_continuous(name="Mean Demand (in kWh)") +
  scale_x_continuous(name="Forecast Uncertainty (in average kWh CRPS)",
                     limits=c(0, 0.1),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))
ggplot3
# ggsave(paste(Sys.Date(),plt3nam,sep="_"),path="./Plots")

plt3rnam <- "rndgrp_01.pdf"
plot3r   <- plot3[plot3$Grouping == 'Random',]
plot3r$uDemand <- round(plot3r$uDemand * 0.4)
dt <- data.table(plot3r)
setkey(dt,uDemand)
plot3rs <- as.data.frame(dt[,mean(CRPS),by=uDemand])
plot3rs$uDemand <- plot3rs$uDemand / 0.4
plot3rs[1,1]=1.25
ggplot3rs <- ggplot(plot3rs, aes(V1,uDemand)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray90"),
                     panel.grid.minor = element_line(colour='gray95'), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),  
        axis.title.x = element_text(color="black",size=18),
        axis.title.y = element_text(color="black",size=18)) +
  labs(x = "Forecast Uncertainty (average kWh CRPS)", y = "Demand (kWh)") +
  scale_x_continuous(expand = c(0, 0),limits = c(0.005,0.055)) + scale_y_continuous(expand = c(0, 0),limits=c(0,26))
ggplot3rs

plt4nam <- "benchKDxAG.pdf"
plot4i  <- readRDS("smuf_temp_compare.rds")
# niceleg <- c(paste('KDE',c(4,24),'h window'),paste('ARMA-GARCH GoF',c(NA,0.01,0.05,0.2,0.5,1.0,2.0)))
plot4   <- plot4i[[1]]
# rownames(plot4) <- niceleg
plot4   <- melt(t(plot4), id=colnames(t(plot4)))
plot4[,1]       <- as.numeric(sub(".*\\.", "", plot4[,1]))
colnames(plot4) <- c('ahead_t','Method','CRPS')
ggplot4 <- ggplot(plot4, aes(ahead_t,CRPS, linetype=Method)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray90"),
                     panel.grid.minor = element_line(colour = "gray95"), axis.line = element_line(colour = "gray60")) +
    theme(text=element_text(family="Times",size=18)) +
    theme(text=element_text(family="Times"),
          axis.text.x = element_text(color="black",size=18),
          axis.text.y = element_text(color="black",size=18),  
          axis.title.x = element_text(color="black",size=18),
          axis.title.y = element_text(color="black",size=18),
          legend.text = element_text(color="black",size=18)) +
    scale_x_continuous(name="Time ahead forecast (h)") +
    scale_y_continuous(limits=c(0.09, 0.121))#,breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
ggplot4
# ggsave(paste(Sys.Date(),plt4nam,sep="_"),path="./Plots")

plt4snam <- "benchKDxAG_simple.pdf"
plot4s <- plot4
levels(plot4s$Method) <- c(levels(plot4s$Method), "ARMA-GARCH(\u03B4=0.01)","ARMA-GARCH(\u03B4=0.05)","ARMA-GARCH(\u03B4=0.20)","ARMA-GARCH(\u03B4=0.00)","KDE")
{plot4s[plot4s=="KD24"]<-"KDE"
plot4s[plot4s=="AG0.00"]<-"ARMA-GARCH(\u03B4=0.00)"
plot4s[plot4s=="AG0.01"]<-"ARMA-GARCH(\u03B4=0.01)"
plot4s[plot4s=="AG0.05"]<-"ARMA-GARCH(\u03B4=0.05)"
plot4s[plot4s=="AG0.20"]<-"ARMA-GARCH(\u03B4=0.20)"}
ggplot4s <- ggplot(plot4s, aes(ahead_t,CRPS, colour=Method, linetype=Method)) + geom_line() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),  
        axis.title.x = element_text(color="black",size=18),
        axis.title.y = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=14),
        legend.text = element_text(color="black",size=14)) +
  theme(legend.position=c(0.8,0.9)) + 
  scale_x_continuous(name="Forecast lead time (h)") +
  scale_y_continuous(name="CRPS (kW)",limits=c(0.09, 0.121)) + 
  scale_colour_grey()
ggplot4s
# ggsave(paste(Sys.Date(),plt4snam,sep="_"),path="./Plots")
expression(beta)
#===========================================
saveRDS(list(plot1,plot2,plot3),file="smuf_aux_plots.rds")
