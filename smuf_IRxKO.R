datKO     <- readRDS("smuf_import-complete.rds")
parKO     <- readRDS("smuf_import-parameter.rds")
datIR     <- readRDS("smuf_import-completeIR.rds")
parIR     <- readRDS("smuf_import-parameterIR.rds")

KOmu <- rowMeans(datKO)
KOsd <- fx_sd_mymat(datKO)
KOmu <- KOmu[-c(710,978)]
KOsd <- KOsd[-c(710,978)]
IRmu <- rowMeans(datIR)
IRsd <- fx_sd_mymat(datIR)

plot(KOmu,KOsd)
plot(IRmu,IRsd)