datKO     <- readRDS("smuf_import-complete.rds")
parKO     <- readRDS("smuf_import-parameter.rds")
datIR     <- readRDS("smuf_import-completeIRhour.rds")
parIR     <- readRDS("smuf_import-parameter.rds")

# only last 6 months of 2010 for IR
datIR3  <- datIR[,2329:6744]
datIR3u <- colMeans(datIR3)
datKO3  <- datKO[,4369:8784]
datKO3u <- colMeans(datKO3)
dates3  <- sort(rep(seq(from=as.Date("2010/7/1"), to=as.Date("2010/12/31"),by='day'),24))
hours3  <- 1:24
dat3f   <- data.frame(date=dates3,hour=hours3,IR=datIR3u,KO=datKO3u)
dat3fl  <- dat3f %>%
  gather()
  

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

