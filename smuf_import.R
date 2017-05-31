#===========================================
# Smart Metering Forecasting using KDE
#
# Author  Estevao "Steve" Alvarenga
#         efsa@bath.edu
# Created in 10/Feb/17
#===========================================
# KSSM0100 importing data
#===========================================

# install.packages("readr")
# install.packages("readxl")

library(readr)

ptm <- proc.time() # Start the clock!

#===========================================
# Parameters
#===========================================
# Period data
s01           = 24          # 1st seazonality cycle
s02           = 168         # 2nd seazonality cycle
s03           = 8760        # 3rd seazonality cycle (not using)
m31           = 31*s01+4    # setting months
m30           = 30*s01+4
m29           = 29*s01+4
m28           = 28*s01+4
y365          = 365*s01     # setting years
y366          = 366*s01
# Selection parameters
sum_of_h      = 1
cus_list      = seq(1,1000)
data_size     = (y366+y365*2-(30+31)*24)/sum_of_h
read_xls_file = T
# Parameter Bundle
importpar     = c(s01,s02,s03,sum_of_h,data_size)

#===========================================
# Read and clean data
#===========================================

# Read data, per excel sheet ===============
if (read_xls_file == TRUE){
  library(readxl)
  KO201201dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 1)
  KO201202dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 2)
  KO201203dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 3)
  KO201204dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 4)
  KO201205dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 5)
  KO201206dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 6)
  KO201207dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 7)
  KO201208dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 8)
  KO201209dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 9)
  KO201210dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 10)
  KO201211dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 11)
  KO201212dom <- read_excel("../data_KO/2012 domestic.xlsx", sheet = 12)

  KO201301dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 1)
  KO201302dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 2)
  KO201303dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 3)
  KO201304dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 4)
  KO201305dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 5)
  KO201306dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 6)
  KO201307dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 7)
  KO201308dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 8)
  KO201309dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 9)
  KO201310dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 10)
  KO201311dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 11)
  KO201312dom <- read_excel("../data_KO/2013 domestic.xlsx", sheet = 12)

  KO201401dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 1)
  KO201402dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 2)
  KO201403dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 3)
  KO201404dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 4)
  KO201405dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 5)
  KO201406dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 6)
  KO201407dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 7)
  KO201408dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 8)
  KO201409dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 9)
  KO201410dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 10)
  # KO201411dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 11)
  # KO201412dom <- read_excel("../data_KO/2014 domestic.xlsx", sheet = 12)
}

# WorkMatrix 01_00 =========================
# Cleared load data for selected customers
wm01_00 = matrix(nrow=(length(cus_list)), ncol=(data_size))
for (j in cus_list){
  if (j == 1){cat("\nLoading working matrix for",length(cus_list),"customers\n")}
  cat(" ",j)
  CUSno       = j
  wv01a       = c(KO201201dom[CUSno,5:m31], KO201202dom[CUSno,5:m29], KO201203dom[CUSno,5:m31], KO201204dom[CUSno,5:m30],
                  KO201205dom[CUSno,5:m31], KO201206dom[CUSno,5:m30], KO201207dom[CUSno,5:m31], KO201208dom[CUSno,5:m31],
                  KO201209dom[CUSno,5:m30], KO201210dom[CUSno,5:m31], KO201211dom[CUSno,5:m30], KO201212dom[CUSno,5:m31])
  wv01b       = c(KO201301dom[CUSno,5:m31], KO201302dom[CUSno,5:m28], KO201303dom[CUSno,5:m31], KO201304dom[CUSno,5:m30],
                  KO201305dom[CUSno,5:m31], KO201306dom[CUSno,5:m30], KO201307dom[CUSno,5:m31], KO201308dom[CUSno,5:m31],
                  KO201309dom[CUSno,5:m30], KO201310dom[CUSno,5:m31], KO201311dom[CUSno,5:m30], KO201312dom[CUSno,5:m31])
  wv01c       = c(KO201401dom[CUSno,5:m31], KO201402dom[CUSno,5:m28], KO201403dom[CUSno,5:m31], KO201404dom[CUSno,5:m30],
                  KO201405dom[CUSno,5:m31], KO201406dom[CUSno,5:m30], KO201407dom[CUSno,5:m31], KO201408dom[CUSno,5:m31],
                  KO201409dom[CUSno,5:m30], KO201410dom[CUSno,5:m31])
  wv01        = c(wv01a,wv01b,wv01c)
  wv02        = as.numeric(wv01)
  for (i in 1:(length(wv02))) {
    if (is.na(wv02[i]) == TRUE) {
      if (i<170){
        wv02[i] = 0
      } else {
        if (is.na(wv02[(i+1)]) == TRUE) {
          wv02[i] = wv02[(i-168)]
        }
        else {wv02[i] = (wv02[i-1]+wv02[i+1])/2}
      }
    }
  }
  wv02_01 = as.numeric(tapply(wv02,(seq_along(wv02)-1) %/% (sum_of_h), sum))
  wm01_00[j,1:(data_size)]  = wv02_01 / sum_of_h
}

saveRDS(wm01_00,   file="smuf_import-complete.rds")
saveRDS(importpar, file="smuf_import-parameter.rds")

print(proc.time() - ptm)        # Stop the clock
