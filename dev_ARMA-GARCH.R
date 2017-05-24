install.packages("rgl")
install.packages("rugarch")
library(rgl)
library(rugarch)
data <- rnorm(1000)
spec <- ugarchspec(variance.model = list(model = "sGARCH", 
                                         garchOrder = c(1, 1), 
                                         submodel = NULL, 
                                         external.regressors = NULL, 
                                         variance.targeting = FALSE), 
                   mean.model     = list(armaOrder = c(1, 1), 
                                         external.regressors = NULL), 
                   distribution.model = "norm", 
                   start.pars = list(), 
                   fixed.pars = list())
spe2 <- ugarchspec(variance.model = list(model = "sGARCH", 
                                         garchOrder = c(2, 2), 
                                         submodel = NULL, 
                                         external.regressors = NULL, 
                                         variance.targeting = FALSE), 
                   mean.model     = list(armaOrder = c(2, 2), 
                                         external.regressors = NULL), 
                   distribution.model = "norm", 
                   start.pars = list(), 
                   fixed.pars = list())

garch <- ugarchfit(spec = spec, data = data, out.sample = 50, solver = "solnp", solver.control = list(trace=0))
garc2 <- ugarchfit(spec = spe2, data = data, out.sample = 50, solver = "solnp", solver.control = list(trace=0))
