#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 13/Sep/17
#===========================================
# KSSM0100 importing data
#===========================================

library(readr)
library(dplyr)
library(tidyr)

#===========================================
# Parameters
#===========================================
# Period data
s01           = 24*2          # 1st seazonality cycle
s02           = 168*2         # 2nd seazonality cycle
s03           = 8760*2        # 3rd seazonality cycle (not using)
sum_of_h      = 2

#===========================================
# Shaping matrix
#===========================================
irishSM <- read_csv("~/Dropbox/PhD student - Steve/Irish Smart Metering Data/1EE customer usage data.csv")
irishSM1 <- irishSM[,2:51]
irishSM1$Date <- unclass(as.Date(irishSM$Date,format='%d/%m/%y'))
irishSM1$Date <- (irishSM1$Date - min(irishSM1$Date)) *s01
irishSM2 <- irishSM1 %>%
  gather(hour,demand,-Date,-Station)
irishSM2$hour <- as.numeric(irishSM2$hour)
irishSM2 <- irishSM2 %>%
  mutate(datehour = Date + hour)
irishSM2$demand <- as.numeric(irishSM2$demand)
irishSM2b <- irishSM2[,c(1,5,4)]
sum(is.na(irishSM2b$demand))
irishSM3 <- irishSM2b %>%
  spread(datehour,demand)

irishSM4 <- as.matrix(irishSM3[,2:13489])
irishSM4 <- unname(irishSM4)
irishSM4 <- t(irishSM4)
length(which(is.na(irishSM4)))
while (length(which(is.na(irishSM4))) > 0) {
  print(length(which(is.na(irishSM4))))
  irishSM4[which(is.na(irishSM4))] = irishSM4[(which(is.na(irishSM4))-1)]
}
length(which(is.na(irishSM4)))
irishSM4 <- t(irishSM4)

irishSM5 <- matrix(ncol=6744,nrow=929)
for (i in 1:nrow(irishSM4)) {
  vect <- as.numeric(tapply(irishSM4[i,],(seq_along(irishSM4[i,])-1) %/% (sum_of_h), sum))
  irishSM5[i,]  = vect / sum_of_h
}

# Parameter bundle
data_size     = ncol(irishSM4)/sum_of_h
importpar     = c(s01,s02,s03,sum_of_h,data_size)

# Saving
saveRDS(irishSM5,  file="smuf_import-completeIRhour.rds")
saveRDS(importpar, file="smuf_import-parameterIR.rds")
