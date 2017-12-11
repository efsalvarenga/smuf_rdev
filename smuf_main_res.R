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
library(tidyr)
library(magrittr)
library(data.table)

fontsize <- 16

setwd("~/GitRepos/smuf_rdev")
source("smuf_main-fxs.R")
{
  # wm01_00       <- readRDS("smuf_import-complete.rds")
  # importpar     <- readRDS("smuf_import-parameter.rds")
  wm01_00       <- readRDS("smuf_import-complete.rds")
  importpar     <- readRDS("smuf_import-parameter.rds")
  implotdt      <- readRDS("smuf_aux_plots.rds")
  s01           <- importpar[1]
  s02           <- importpar[2]
  s03           <- importpar[3]
  sum_of_h      <- importpar[4]
  data_size     <- importpar[5]
}

# compare_demands
cus_no  <- 900
ds_ini  <- 2001
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
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
ggplot1b <- ggplot(plot1, aes(Time,Cus4)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
ggplot1c <- ggplot(plot1, aes(Time,Cus5)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
ggplot1d <- ggplot(plot1, aes(Time,Cus6)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
multiplot(ggplot1a, ggplot1d, ggplot1c, ggplot1b, cols=2)

# compare_demands_grouping
cus_no  <- 900
ds_ini  <- 2001
ds_len  <- s02*4
plot1B  <- as.data.frame(t(wm01_00[1:cus_no,ds_ini:(ds_ini+ds_len-1)]))
plot1B  <- cbind(seq(1,ds_len),plot1B)
colnames(plot1B) <- c("Time",paste("Cus",1:cus_no,sep=""))
plot1C1 <- plot1B[,c(1,3)]
plot1C2 <- rowSums(plot1B[,2:11])/10
plot1C3 <- rowSums(plot1B[,2:101])/100
plot1C4 <- rowSums(plot1B[,2:201])/200
plot1C  <- cbind(plot1C1,plot1C2,plot1C3,plot1C4)
colnames(plot1C) <- c("Time",'a1cus','a10cus','a100cus','a200cus')
ggplot1Ca <- ggplot(plot1C, aes(Time,a1cus)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
ggplot1Cb <- ggplot(plot1C, aes(Time,a10cus)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
ggplot1Cc <- ggplot(plot1C, aes(Time,a100cus)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
ggplot1Cd <- ggplot(plot1C, aes(Time,a200cus)) +
  geom_line(color="black") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  labs(x = "Hours", y = "Demand (kWh)") +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize*1.2),
        axis.text.y = element_text(color="black",size=fontsize*1.2),  
        axis.title.x = element_text(color="black",size=fontsize*1.2),
        axis.title.y = element_text(color="black",size=fontsize*1.2))
multiplot(ggplot1Ca, ggplot1Cc, ggplot1Cb, ggplot1Cd, cols=2)

# compare_densities.pdf
# plot2   <- implotdt[[2]]
{cus_nos <- c(8,9,12,14) # for KO
# cus_nos <- c(2,5,12,16)#seq(11,15)) # for IR
ds_init <- 1993
sl_win  <- seq(ds_init,(ds_init+s01-1))
plot2   <- as.numeric(wm01_00[cus_nos,sl_win])
plot2   <- as.data.frame(plot2)
plot2   <- cbind(plot2,rep(paste("Cus",cus_nos,sep=""),s01))
colnames(plot2) <- c("Demand","Customer")}
ggplot2 <- ggplot(plot2, aes(Demand, linetype=Customer)) +
            theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
            geom_density(alpha = 0.1) +
            theme(text=element_text(family="Times"),
                  axis.text.x = element_text(color="black",size=fontsize),
                  axis.text.y = element_text(color="black",size=fontsize),  
                  axis.title.x = element_text(color="black",size=fontsize),
                  axis.title.y = element_text(color="black",size=fontsize)) +
            theme(legend.position="none") +
            labs(x = "Demand (kWh)", y = "Density") +
            scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
ggplot2

# plt3nam <- "optgrp_01.pdf"
# # plot3i  <- readRDS("smuf_run_0624_defheurd.rds")
plot3i  <- readRDS("smuf_run_0704_defheur_med01.rds")
myleg   <- c("Random","SDKD","SDAG","CVKD","CVAG")
# fx_plt_rnd_vs_opt(plot3i[[2]][[length(plot3i[[2]])]][[2]],c(0,0.05),c(0,27),myleg,"CRPS")
plot3   <- cbind(as.data.frame(plot3i[[2]][[1]][[2]][[2]]),myleg[1])
colnames(plot3) <- c("CRPS","uDemand","Grouping")
# niceleg <- c("Optimal (st. dev.)","","Optimal (cross-val.)")
# for (i in c(4,6)) {
#   temp  <- cbind(as.data.frame(plot3i[[2]][[1]][[2]][[i]]),niceleg[(i-3)])
#   colnames(temp) <- c("CRPS","uDemand","Grouping")
#   plot3 <- rbind(plot3,temp)
# }
# ggplot3 <- ggplot(plot3, aes(CRPS,uDemand, color=Grouping)) + geom_point() +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray90"),
#                      panel.grid.minor = element_line(colour = "gray95"), axis.line = element_line(colour = "gray60")) +
#   scale_color_manual(values=c("gray80", "dodgerblue3", "firebrick")) +
#   theme(text=element_text(family="Times",size=fontsize)) +
#   scale_y_continuous(name="Mean Demand (in kWh)") +
#   scale_x_continuous(name="Forecast Uncertainty (in average kWh CRPS)",
#                      limits=c(0, 0.1),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
#   theme(legend.position=c(0.8,0.7))
# ggplot3
# ggsave(paste(Sys.Date(),plt3nam,sep="_"),path="./Plots")

plt3rnam <- "rndgrp_01.pdf"
plot3r   <- plot3[plot3$Grouping == 'Random',]
plot3r$uDemand <- round(plot3r$uDemand * 0.4)
dt <- data.table(plot3r)
setkey(dt,uDemand)
plot3rs <- as.data.frame(dt[,mean(CRPS),by=uDemand])
plot3rs$uDemand <- plot3rs$uDemand / 0.4
plot3rs[1,1]=1.25
ggplot3rs <- ggplot(plot3rs, aes(V1,uDemand)) + geom_point(col='gray30') +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize)) +
  labs(x = "CRPS (in kWh)", y = "Mean Demand (kWh)") +
  scale_x_continuous(expand = c(0, 0),limits = c(0.005,0.055)) + scale_y_continuous(expand = c(0, 0),limits=c(0,26))
ggplot3rs

# plt4nam <- "benchKDxAG.pdf"
# plot4i  <- readRDS("smuf_temp_compare.rds")
# # niceleg <- c(paste('KDE',c(4,24),'h window'),paste('ARMA-GARCH GoF',c(NA,0.01,0.05,0.2,0.5,1.0,2.0)))
# plot4   <- plot4i[[1]]
# # rownames(plot4) <- niceleg
# plot4   <- melt(t(plot4), id=colnames(t(plot4)))
# plot4[,1]       <- as.numeric(sub(".*\\.", "", plot4[,1]))
# colnames(plot4) <- c('ahead_t','Method','CRPS')
# ggplot4 <- ggplot(plot4, aes(ahead_t,CRPS, linetype=Method)) + geom_line() +
#     theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "gray90"),
#                      panel.grid.minor = element_line(colour = "gray95"), axis.line = element_line(colour = "gray60")) +
#     theme(text=element_text(family="Times",size=fontsize)) +
#     theme(text=element_text(family="Times"),
#           axis.text.x = element_text(color="black",size=fontsize),
#           axis.text.y = element_text(color="black",size=fontsize),  
#           axis.title.x = element_text(color="black",size=fontsize),
#           axis.title.y = element_text(color="black",size=fontsize),
#           legend.text = element_text(color="black",size=fontsize)) +
#     scale_x_continuous(name="Time ahead forecast (h)") +
#     scale_y_continuous(limits=c(0.09, 0.121))#,breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
# ggplot4
# # ggsave(paste(Sys.Date(),plt4nam,sep="_"),path="./Plots")

# benchKDxAG_simple
plot4s <- readRDS('smuf_compare_0831_KD-AG_large_mult-gofmin_cleaned.rds')
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
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_text(color="black",size=fontsize),
        legend.text = element_text(color="black",size=fontsize)) +
  theme(legend.position=c(0.8,0.85)) + 
  scale_x_continuous(name="Forecast lead time (h)") +
  scale_y_continuous(name="CRPS (in kWh)",limits=c(0.09, 0.121)) + 
  scale_colour_grey()
ggplot4s

# # seaf24
# plot5A    <- readRDS("smuf_runf_0919_KO_seaf10-summary.rds")
# correc <- (max(plot5A[[2]][,2])-min(plot5A[[2]][,2]))/16/2
# corvec <- runif(1600,-correc,correc)
# plot5A[[2]][,2] <- plot5A[[2]][,2] + corvec
# myleg     <- c("Random","Standard Deviation","Seasonal Signal","Seas. & Remainder Signal")
# plot5A[[2]] <- as.data.frame(cbind(plot5A[[2]],myleg[1]))
# colnames(plot5A[[2]]) <- c("CRPS","uDemand","Grouping")
# plot5A[[3]] <- as.data.frame(cbind(plot5A[[3]],myleg[2]))
# plot5A[[4]] <- as.data.frame(cbind(plot5A[[4]],myleg[3]))
# plot5A[[5]] <- as.data.frame(cbind(plot5A[[5]],myleg[4]))
# colnames(plot5A[[2]]) <- c("CRPS","uDemand","Grouping")
# colnames(plot5A[[3]]) <- c("CRPS","uDemand","Grouping")
# colnames(plot5A[[4]]) <- c("CRPS","uDemand","Grouping")
# colnames(plot5A[[5]]) <- c("CRPS","uDemand","Grouping")
# plot5A    <- rbind(plot5A[[2]],plot5A[[3]],plot5A[[4]],plot5A[[5]],make.row.names = FALSE)
# plot5A$CRPS <- as.numeric(as.character(plot5A$CRPS))
# plot5A$uDemand <- as.numeric(as.character(plot5A$uDemand))
# 
# ggplot5A  <- ggplot(plot5A, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
#   theme(text=element_text(family="Times"),
#         axis.text.x = element_text(color="black",size=fontsize),
#         axis.text.y = element_text(color="black",size=fontsize),  
#         axis.title.x = element_text(color="black",size=fontsize),
#         axis.title.y = element_text(color="black",size=fontsize),
#         legend.title = element_blank(),
#         legend.text = element_text(color="black",size=fontsize)) +
#   scale_color_manual(values=c("gray80", "black", "black", "black")) +
#   scale_y_continuous(name="Mean Demand (in kWh)") +
#   scale_x_continuous(name="CRPS (in kWh)",
#                      limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
#   theme(legend.position=c(0.8,0.7))
# ggplot5A
# 
# plot5B    <- readRDS("smuf_runf_0919_KO_base_cont_real.rds")
# plot5Brnd <- plot5B[[2]][[1]][[2]][[2]]
# plot5Bsd  <- plot5B[[2]][[1]][[2]][[3]]
# plot5Bcv  <- plot5B[[2]][[1]][[2]][[4]]
# for (i in c(2:10)) {
#   plot5Brnd <- plot5Brnd + plot5B[[2]][[i]][[2]][[2]]
#   plot5Bsd  <- plot5Bsd  + plot5B[[2]][[i]][[2]][[3]]
#   plot5Bcv  <- plot5Bcv  + plot5B[[2]][[i]][[2]][[4]]
# }
# plot5Brnd <- plot5Brnd/10
# plot5Bsd  <- plot5Bsd/10
# plot5Bcv  <- plot5Bcv/10
# correc <- (max(plot5Brnd[,2])-min(plot5Brnd[,2]))/16/2
# corvec <- runif(1600,-correc,correc)
# plot5Brnd[,2] <- plot5Brnd[,2] + corvec
# myleg     <- c("Random","Standard Deviation","Forecast Validated")
# plot5Brnd <- as.data.frame(cbind(plot5Brnd,myleg[1]))
# plot5Bsd  <- as.data.frame(cbind(plot5Bsd,myleg[2]))
# plot5Bcv  <- as.data.frame(cbind(plot5Bcv,myleg[3]))
# colnames(plot5Brnd) <- c("CRPS","uDemand","Grouping")
# colnames(plot5Bsd)  <- c("CRPS","uDemand","Grouping")
# colnames(plot5Bcv)  <- c("CRPS","uDemand","Grouping")
# plot5B      <- rbind(plot5Brnd,plot5Bsd,plot5Bcv,make.row.names = FALSE)
# plot5B$CRPS <- as.numeric(as.character(plot5B$CRPS))
# plot5B$uDemand <- as.numeric(as.character(plot5B$uDemand))
# 
# ggplot5B  <- ggplot(plot5B, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
#   theme(text=element_text(family="Times"),
#         axis.text.x = element_text(color="black",size=fontsize),
#         axis.text.y = element_text(color="black",size=fontsize),  
#         axis.title.x = element_text(color="black",size=fontsize),
#         axis.title.y = element_text(color="black",size=fontsize),
#         legend.title = element_blank(),
#         legend.text = element_text(color="black",size=fontsize)) +
#   scale_color_manual(values=c("gray80", "black", "black")) +
#   scale_y_continuous(name="Mean Demand (in kWh)") +
#   scale_x_continuous(name="CRPS (in kWh)",
#                      limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
#   theme(legend.position=c(0.8,0.7))
# ggplot5B
# 
# plot5AB <- plot5A %>%
#   filter(Grouping %in% c("Seasonal Signal","Seas. & Remainder Signal"))
# plot5AB <- rbind(plot5B,plot5AB)
# ggplot5AB  <- ggplot(plot5AB, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
#   theme(text=element_text(family="Times"),
#         axis.text.x = element_text(color="black",size=fontsize),
#         axis.text.y = element_text(color="black",size=fontsize),  
#         axis.title.x = element_text(color="black",size=fontsize),
#         axis.title.y = element_text(color="black",size=fontsize),
#         legend.title = element_blank(),
#         legend.text = element_text(color="black",size=fontsize)) +
#   scale_color_manual(values=c("gray80", rep("black", 4))) +
#   scale_y_continuous(name="Mean Demand (in kWh)") +
#   scale_x_continuous(name="CRPS (kW)",
#                      limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
#   theme(legend.position=c(0.8,0.7))
# ggplot5AB

# seaf04
plot6A    <- readRDS("smuf_runf_0919_KO_seaf10-4h_summary.rds")
correc <- (max(plot6A[[2]][,2])-min(plot6A[[2]][,2]))/16/2
corvec <- runif(1600,-correc,correc)
plot6A[[2]][,2] <- plot6A[[2]][,2] + corvec
myleg     <- c("Random","Standard Deviation","Seasonal Signal","Seas. & Remainder Signal")
plot6A[[2]] <- as.data.frame(cbind(plot6A[[2]],myleg[1]))
colnames(plot6A[[2]]) <- c("CRPS","uDemand","Grouping")
plot6A[[3]] <- as.data.frame(cbind(plot6A[[3]],myleg[2]))
plot6A[[4]] <- as.data.frame(cbind(plot6A[[4]],myleg[3]))
plot6A[[5]] <- as.data.frame(cbind(plot6A[[5]],myleg[4]))
colnames(plot6A[[2]]) <- c("CRPS","uDemand","Grouping")
colnames(plot6A[[3]]) <- c("CRPS","uDemand","Grouping")
colnames(plot6A[[4]]) <- c("CRPS","uDemand","Grouping")
colnames(plot6A[[5]]) <- c("CRPS","uDemand","Grouping")
plot6A    <- rbind(plot6A[[2]],plot6A[[3]],plot6A[[4]],plot6A[[5]],make.row.names = FALSE)
plot6A$CRPS <- as.numeric(as.character(plot6A$CRPS))
plot6A$uDemand <- as.numeric(as.character(plot6A$uDemand))

ggplot6A  <- ggplot(plot6A, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=fontsize)) +
  scale_color_manual(values=c("gray80", "black", "black", "black")) +
  scale_y_continuous(name="Mean Demand (in kWh)") +
  scale_x_continuous(name="CRPS (in kWh)",
                     limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))
ggplot6A

plot6B    <- readRDS("smuf_runf_0919_KO_base-4h_summary.rds")
correc <- (max(plot6B[[2]][,2])-min(plot6B[[2]][,2]))/16/2
corvec <- runif(1600,-correc,correc)
plot6B[[2]][,2] <- plot6B[[2]][,2] + corvec
myleg     <- c("Random","Standard Deviation","Forecast Validated")
plot6B[[2]] <- as.data.frame(cbind(plot6B[[2]],myleg[1]))
colnames(plot6B[[2]]) <- c("CRPS","uDemand","Grouping")
plot6B[[3]] <- as.data.frame(cbind(plot6B[[3]],myleg[2]))
plot6B[[4]] <- as.data.frame(cbind(plot6B[[4]],myleg[3]))
colnames(plot6B[[2]]) <- c("CRPS","uDemand","Grouping")
colnames(plot6B[[3]]) <- c("CRPS","uDemand","Grouping")
colnames(plot6B[[4]]) <- c("CRPS","uDemand","Grouping")
plot6B    <- rbind(plot6B[[2]],plot6B[[3]],plot6B[[4]],make.row.names = FALSE)
plot6B$CRPS <- as.numeric(as.character(plot6B$CRPS))
plot6B$uDemand <- as.numeric(as.character(plot6B$uDemand))

plot6AB <- plot6A %>%
  filter(Grouping %in% c("Seasonal Signal","Seas. & Remainder Signal"))
plot6AB <- rbind(plot6B,plot6AB)
ggplot6AB  <- ggplot(plot6AB, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=fontsize)) +
  scale_color_manual(values=c("gray80", rep("black", 4))) +
  scale_y_continuous(name="Mean Demand (in kWh)") +
  scale_x_continuous(name="CRPS (kW)",
                     limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))
ggplot6AB

# all12
{plot7    <- readRDS("smuf_runf_1121_KO_all-12h_ws04_summary.rds")
correc <- (max(plot7[[2]][,2])-min(plot7[[2]][,2]))/16/2
corvec <- runif(1600,-correc,correc)
plot7[[2]][,2] <- plot7[[2]][,2] + corvec
myleg     <- c("Random","Seasonal Residual","Forecast Validated","Seasonal Similarity","Seas. & Remainder")
plot7[[2]] <- as.data.frame(cbind(plot7[[2]],myleg[1]))
colnames(plot7[[2]]) <- c("CRPS","uDemand","Grouping")
plot7[[3]] <- as.data.frame(cbind(plot7[[3]],myleg[2]))
plot7[[4]] <- as.data.frame(cbind(plot7[[4]],myleg[3]))
plot7[[5]] <- as.data.frame(cbind(plot7[[5]],myleg[4]))
plot7[[6]] <- as.data.frame(cbind(plot7[[6]],myleg[5]))
colnames(plot7[[2]]) <- c("CRPS","uDemand","Grouping")
colnames(plot7[[3]]) <- c("CRPS","uDemand","Grouping")
colnames(plot7[[4]]) <- c("CRPS","uDemand","Grouping")
colnames(plot7[[5]]) <- c("CRPS","uDemand","Grouping")
colnames(plot7[[6]]) <- c("CRPS","uDemand","Grouping")
plot7    <- rbind(plot7[[2]],plot7[[3]],plot7[[4]],plot7[[5]],plot7[[6]],make.row.names = FALSE)
plot7$CRPS <- as.numeric(as.character(plot7$CRPS))
plot7$uDemand <- as.numeric(as.character(plot7$uDemand))
ggplot7  <- ggplot(plot7, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=fontsize)) +
  scale_color_manual(values=c("gray80", rep("black", 4))) +
  scale_y_continuous(name="Mean Demand (in kWh)") +
  scale_x_continuous(name="CRPS (kW)",
                     limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))}
ggplot7

# all04
{plot8    <- readRDS("smuf_runf_1121_KO_all-04h_ws04_summary.rds")
correc <- (max(plot8[[2]][,2])-min(plot8[[2]][,2]))/16/2
corvec <- runif(1600,-correc,correc)
plot8[[2]][,2] <- plot8[[2]][,2] + corvec
myleg     <- c("Random","Seasonal Residual","Forecast Validated","Seasonal Similarity","Seas. & Remainder")
plot8[[2]] <- as.data.frame(cbind(plot8[[2]],myleg[1]))
colnames(plot8[[2]]) <- c("CRPS","uDemand","Grouping")
plot8[[3]] <- as.data.frame(cbind(plot8[[3]],myleg[2]))
plot8[[4]] <- as.data.frame(cbind(plot8[[4]],myleg[3]))
plot8[[5]] <- as.data.frame(cbind(plot8[[5]],myleg[4]))
plot8[[6]] <- as.data.frame(cbind(plot8[[6]],myleg[5]))
colnames(plot8[[2]]) <- c("CRPS","uDemand","Grouping")
colnames(plot8[[3]]) <- c("CRPS","uDemand","Grouping")
colnames(plot8[[4]]) <- c("CRPS","uDemand","Grouping")
colnames(plot8[[5]]) <- c("CRPS","uDemand","Grouping")
colnames(plot8[[6]]) <- c("CRPS","uDemand","Grouping")
plot8    <- rbind(plot8[[2]],plot8[[3]],plot8[[4]],plot8[[5]],plot8[[6]],make.row.names = FALSE)
plot8$CRPS <- as.numeric(as.character(plot8$CRPS))
plot8$uDemand <- as.numeric(as.character(plot8$uDemand))
ggplot8  <- ggplot(plot8, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=fontsize)) +
  scale_color_manual(values=c("gray80", rep("black", 4))) +
  scale_y_continuous(name="Mean Demand (in kWh)") +
  scale_x_continuous(name="CRPS (kW)",
                     limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))}
ggplot8

# rlx04
{plot9    <- readRDS("smuf_runf_1121_KO_all-04h_ws04_rlx_summary.rds")
correc <- (max(plot9[[2]][,2])-min(plot9[[2]][,2]))/16/2
corvec <- runif(1600,-correc,correc)
plot9[[2]][,2] <- plot9[[2]][,2] + corvec
myleg     <- c("Random","Seasonal Residual","SD Binary Relaxed")
plot9[[2]] <- as.data.frame(cbind(plot9[[2]],myleg[1]))
colnames(plot9[[2]]) <- c("CRPS","uDemand","Grouping")
plot9[[3]] <- as.data.frame(cbind(plot9[[3]],myleg[2]))
plot9[[4]] <- as.data.frame(cbind(plot9[[4]],myleg[3]))
colnames(plot9[[2]]) <- c("CRPS","uDemand","Grouping")
colnames(plot9[[3]]) <- c("CRPS","uDemand","Grouping")
colnames(plot9[[4]]) <- c("CRPS","uDemand","Grouping")
plot9    <- rbind(plot9[[2]],plot9[[3]],plot9[[4]],make.row.names = FALSE)
plot9$CRPS <- as.numeric(as.character(plot9$CRPS))
plot9$uDemand <- as.numeric(as.character(plot9$uDemand))
ggplot9  <- ggplot(plot9, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=fontsize)) +
  scale_color_manual(values=c("gray80", rep("black", 4))) +
  scale_y_continuous(name="Mean Demand (in kWh)") +
  scale_x_continuous(name="CRPS (kW)",
                     limits=c(0, 0.08),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
  theme(legend.position=c(0.8,0.7))}
ggplot9

# irish
{plot10    <- readRDS("smuf_runf_1121_IR_sd-04h_ws04_summary.rds")
  correc <- (max(plot10[[2]][,2])-min(plot10[[2]][,2]))/16/2
  corvec <- runif(1600,-correc,correc)
  plot10[[2]][,2] <- plot10[[2]][,2] + corvec
  myleg     <- c("Random","Seasonal Residual")
  plot10[[2]] <- as.data.frame(cbind(plot10[[2]],myleg[1]))
  colnames(plot10[[2]]) <- c("CRPS","uDemand","Grouping")
  plot10[[3]] <- as.data.frame(cbind(plot10[[3]],myleg[2]))
  colnames(plot10[[2]]) <- c("CRPS","uDemand","Grouping")
  colnames(plot10[[3]]) <- c("CRPS","uDemand","Grouping")
  plot10    <- rbind(plot10[[2]],plot10[[3]],make.row.names = FALSE)
  plot10$CRPS <- as.numeric(as.character(plot10$CRPS))
  plot10$uDemand <- as.numeric(as.character(plot10$uDemand))
  ggplot10  <- ggplot(plot10, aes(CRPS,uDemand, color=Grouping, shape=Grouping)) + geom_point() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
    theme(text=element_text(family="Times"),
          axis.text.x = element_text(color="black",size=fontsize),
          axis.text.y = element_text(color="black",size=fontsize),  
          axis.title.x = element_text(color="black",size=fontsize),
          axis.title.y = element_text(color="black",size=fontsize),
          legend.title = element_blank(),
          legend.text = element_text(color="black",size=fontsize)) +
    scale_color_manual(values=c("gray80", rep("black", 4))) +
    scale_y_continuous(name="Mean Demand (in kWh)") +
    scale_x_continuous(name="CRPS (kW)",
                       limits=c(0, 0.04),breaks=seq(0,0.1,0.02)) +#,expand=c(0,0)) +
    theme(legend.position=c(0.8,0.7))}
ggplot10

plot7 %>% #12h
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS))

plot8 %>% #04h
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS))

plot9 %>% #04h-rlx
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS))

plot10 %>% #irish
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS))

# random lines
plot11     <- readRDS('smuf_runf_0919_KO_randomlines.rds')
plot11     <- t(plot11)
rnd.names <- c(1,2,3,4,5,10,20,30,40,50,100,150,200)
plot11     <- cbind(1:72,plot11)
colnames(plot11) <- c('ahead_t',rnd.names)
plot11     <- as.data.frame(plot11)
plot11     <- plot11 %>% gather(Aggregation,CRPS,-ahead_t)
plot11s    <- plot11 %>% filter (Aggregation %in% c(1,10,100,200))

ggplot11 <- ggplot(plot11s, aes(ahead_t,CRPS, colour=Aggregation, linetype=Aggregation)) + geom_line() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray60")) +
  theme(text=element_text(family="Times"),
        axis.text.x = element_text(color="black",size=fontsize),
        axis.text.y = element_text(color="black",size=fontsize),  
        axis.title.x = element_text(color="black",size=fontsize),
        axis.title.y = element_text(color="black",size=fontsize),
        legend.title = element_text(color="black",size=fontsize),
        legend.text = element_text(color="black",size=fontsize)) +
  theme(legend.position=c(0.92,0.65)) + 
  scale_x_continuous(name="Forecast lead time (h)") +
  scale_y_continuous(name="CRPS (kW)",limits=c(0, 0.115)) + 
  scale_colour_grey()
ggplot11

# compare cv1d or cv7d
table01 <- readRDS('smuf_run_0801_defheur_cv1d_summary.rds')
tb01leg <- c("OptCV7dKD","OptCV1dKD","OptCV7dAG","OptCV1dAG")
table01[[4]] <- as.data.frame(cbind(table01[[4]],tb01leg[1]))
table01[[5]] <- as.data.frame(cbind(table01[[5]],tb01leg[2]))
table01[[6]] <- as.data.frame(cbind(table01[[6]],tb01leg[3]))
table01[[7]] <- as.data.frame(cbind(table01[[7]],tb01leg[4]))
table01b     <- rbind(table01[[4]],table01[[5]],table01[[6]],table01[[7]])
colnames(table01b) <- c("CRPS","uDemand","Grouping")
table01b$CRPS      <- as.numeric(as.character(table01b$CRPS))
table01b$uDemand   <- as.numeric(as.character(table01b$uDemand))
table01b %>%
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS))

# compare seafs
table02 <- readRDS('smuf_run_0720_defheur_clulight_int01_summary.rds')
tb02leg <- c("OptSEAF_nf","OptSEAF_sf","OptSEAF_na","OptSEAF_sa")
table02[[4]] <- as.data.frame(cbind(table02[[4]],tb02leg[1]))
table02[[5]] <- as.data.frame(cbind(table02[[5]],tb02leg[2]))
table02[[6]] <- as.data.frame(cbind(table02[[6]],tb02leg[3]))
table02[[7]] <- as.data.frame(cbind(table02[[7]],tb02leg[4]))
table02b     <- rbind(table02[[4]],table02[[5]],table02[[6]],table02[[7]])
colnames(table02b) <- c("CRPS","uDemand","Grouping")
table02b$CRPS      <- as.numeric(as.character(table02b$CRPS))
table02b$uDemand   <- as.numeric(as.character(table02b$uDemand))
table02b %>%
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS,na.rm=T))

table03 <- readRDS('smuf_run_0720_defheur_clulight_int02_summary.rds')
tb03leg <- c("OptSEAF_na","OptSEAF_na02","OptSEAF_na05","OptSEAF_na08")
table03[[4]] <- as.data.frame(cbind(table03[[4]],tb03leg[1]))
table03[[5]] <- as.data.frame(cbind(table03[[5]],tb03leg[2]))
table03[[6]] <- as.data.frame(cbind(table03[[6]],tb03leg[3]))
table03[[7]] <- as.data.frame(cbind(table03[[7]],tb03leg[4]))
table03b     <- rbind(table03[[4]],table03[[5]],table03[[6]],table03[[7]])
colnames(table03b) <- c("CRPS","uDemand","Grouping")
table03b$CRPS      <- as.numeric(as.character(table03b$CRPS))
table03b$uDemand   <- as.numeric(as.character(table03b$uDemand))
table03b %>%
  group_by(Grouping) %>%
  summarise(uCRPS = mean(CRPS,na.rm=T))

