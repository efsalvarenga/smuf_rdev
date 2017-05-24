library(quantmod)
library(lattice)
library(timeSeries)
library(rugarch)
library(readr)

X_GSPC    <- read_csv("~/Desktop/^GSPC.csv")
spReturns <- diff(log(X_GSPC$Close))
spReturns <- spReturns[!is.na(spReturns)]

windowLength <- 500
foreLength   <- length(spReturns) - windowLength
forecasts    <- vector(mode="character", length=foreLength)

# here starts the loop for rolling forecasts [not implemented loop]
spReturnsOffset = spReturns[(1+d):(windowLength+d)]
print(d)

# finding ARMA spec (p,q -> 0 : 5)
final.aic <- Inf
final.order <- c(0,0,0)
for (p in 0:5) for (q in 0:5) {
  if ( p == 0 && q == 0) {
    next
  }
  arimaFit = tryCatch(arima(spReturnsOffset, order=c(p, 0, q)),
                      error=function( err ) FALSE,
                      warning=function( err ) FALSE )
  if( !is.logical( arimaFit ) ) {
    current.aic <- AIC(arimaFit)
    if (current.aic < final.aic) {
      final.aic <- current.aic
      final.order <- c(p, 0, q)
      final.arima <- arima(spReturnsOffset, order=final.order)
    }
  } else {
    next
  }
}

# setting GARCH spec (1,1)
spec = ugarchspec(variance.model=list(garchOrder=c(1,1)),
                  mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
                  distribution.model="sged")

# fitting ARMA-GARCH parameters
fit = tryCatch(ugarchfit(spec, spReturnsOffset, solver = 'hybrid'),
               error=function(e) e, warning=function(w) w)

# fore = ugarchforecast(fit, n.ahead=3)                       # for point forecast
# dist = ugarchdistribution(fit, n.sim = 100, m.sim = 100)    # for parameters density distr
sim1 = ugarchsim(fit, n.sim = 12, m.sim = 100)                # density simulation forecast

# based on: https://www.quantstart.com/articles/ARIMA-GARCH-Trading-Strategy-on-the-SP500-Stock-Market-Index-Using-R