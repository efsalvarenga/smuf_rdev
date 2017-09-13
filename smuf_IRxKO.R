datKO     <- readRDS("smuf_import-complete.rds")
parKO     <- readRDS("smuf_import-parameter.rds")
datIR     <- readRDS("smuf_import-completeIR.rds")
parIR     <- readRDS("smuf_import-parameterIR.rds")

rowMeans(datKO)
rowMeans(datIR)
sum(is.na(datIR[37,]))
