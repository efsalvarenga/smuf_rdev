datKO     <- readRDS("smuf_import-complete.rds")
parKO     <- readRDS("smuf_import-parameter.rds")
datIR     <- readRDS("smuf_import-completeIR.rds")
parIR     <- readRDS("smuf_import-parameterIR.rds")

datIR2 <- foreach (i = 1:(nrow(datIR)),.combine=c("rbind")) %dopar% {
  unname(tapply(datIR[i,], (seq_along(datIR[i,])-1) %/% 2, sum))
}
datIR2 <- datIR2/2

dim(datKO)
dim(datIR)
dim(datIR2)

KOmu  <- rowMeans(datKO)
KOsd  <- fx_sd_mymat(datKO)
KOmu  <- KOmu[-c(710,979)]
KOsd  <- KOsd[-c(710,979)]
IRmu  <- rowMeans(datIR)
IRsd  <- fx_sd_mymat(datIR)
IRmu  <- IRmu[-c(7)]
IRsd  <- IRsd[-c(7)]
IR2mu <- rowMeans(datIR2)
IR2sd <- fx_sd_mymat(datIR2)
IR2mu  <- IR2mu[-c(7)]
IR2sd  <- IR2sd[-c(7)]

plot(KOmu,KOsd)
plot(IRmu,IRsd)
plot(IR2mu,IR2sd)

