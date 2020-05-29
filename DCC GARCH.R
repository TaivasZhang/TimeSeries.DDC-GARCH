library(fpp)
library(TSA)
library(xts)
library(FinTS)
library(rugarch)
library(rmgarch)
library(quantmod)
library(DataCombine)
library(PerformanceAnalytics)

crude_oil = "Oil Price weekly.csv"

data <- read.csv(crude_oil)
names(data)

data$Date = as.Date(data$ï..Date, "%Y-%m-%d")

par(mfrow=c(2,2))

WTI <- xts(data$WTI, data$Date)
ts.plot(WTI)

BRENT <- xts(data$BRENT, data$Date)
ts.plot(BRENT)

DUBAI <- xts(data$DUBAI, data$Date)
ts.plot(DUBAI)

ESPO <- xts(data$ESPO, data$Date)
ts.plot(ESPO)

MARKET <- xts(cbind(WTI, BRENT, DUBAI, ESPO), index(WTI))
par(mfrow=c(1,1))
plot(
  MARKET, 
  main="Spot Prices of the Major Oil Markets",
  major.ticks="years", minor.ticks=NULL, 
  grid.ticks.on="years",
  col=c("#DE7A22", "#F4CC70", "#20948B", "#6AB187")
)


#Calcuate the log-return of all major markets
market <- log(MARKET)
#When we take the difference, the first element of the original table will be zero.
logreturn <- diff(market)[-1, ]
length(logreturn$WTI)
return.train = logreturn[c(1:360),]

#Test Stationarity
par(mfrow=c(4,1))

# WTI
ts.plot(return.train$WTI,
        xlab="04/19/2013 - 03/06/2020",
        ylab="WTI", col="#CE5A57")
adf.test(return.train$WTI)
kpss.test(return.train$WTI)

# BRENT
ts.plot(return.train$BRENT,
        xlab="04/19/2013 - 03/06/2020",
        ylab="BRENT",col="#F4CC70")
adf.test(return.train$BRENT)
kpss.test(return.train$BRENT)

# Dubai
ts.plot(return.train$DUBAI,
        xlab="04/19/2013 - 03/06/2020",
        ylab="DUBAI",col="#444C5C")
adf.test(return.train$DUBAI)
kpss.test(return.train$DUBAI)

# ESPO
ts.plot(return.train$ESPO,
        xlab="04/19/2013 - 03/06/2020",
        ylab="ESPO", col="#6AB187")
adf.test(return.train$ESPO)
kpss.test(return.train$ESPO)


#Ljung-Box test and ARMA models

for(i in 1:4) print(Box.test(return.train$WTI, lag = 5*i, type =  "Ljung-Box"))
for(i in 1:4) print(Box.test(return.train$BRENT, lag = 5*i, type =  "Ljung-Box"))
for(i in 1:4)print(Box.test(return.train$DUBAI, lag = 5*i, type =  "Ljung-Box"))
for(i in 1:4) print(Box.test(return.train$ESPO, lag = 5*i, type =  "Ljung-Box"))

#Mean Process Filtering using ARMA
wti_arma <- auto.arima(return.train$WTI, stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##wti ARMA(2,1)
brent_arma <- auto.arima(return.train$BRENT, stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##brent ARMA(2,1)
dubai_arma <- auto.arima(return.train$DUBAI, stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##dubai ARMA(0,1) or ARMA(1,1)
espo_arma <- auto.arima(return.train$ESPO, stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##espo ARMA(1,1) or ARMA(0,1)

#Check the distributions of the residuals
# WTI
checkresiduals(wti_arma)
shapiro.test(wti_arma$residuals)
##The residuals are not normally distributed.
# BRENT
checkresiduals(brent_arma)
shapiro.test(brent_arma$residuals)
##The residuals are not normally distributed.
# DUBAI
checkresiduals(dubai_arma)
shapiro.test(dubai_arma$residuals)
##The residuals are not normally distributed.
# ESPO
checkresiduals(espo_arma)
shapiro.test(espo_arma$residuals)
##The residuals are not normally distributed.

for(i in 1:4) print(Box.test(wti_arma$residuals, lag = 5*i, type =  "Ljung-Box"))
for(i in 1:4) print(Box.test(brent_arma$residuals, lag = 5*i, type =  "Ljung-Box"))
for(i in 1:4) print(Box.test(dubai_arma$residuals, lag = 5*i, type =  "Ljung-Box"))
for(i in 1:4) print(Box.test(espo_arma$residuals, lag = 5*i, type =  "Ljung-Box"))

r_wti <- wti_arma$residuals
r_brent <- brent_arma$residuals
r_dubai <- dubai_arma$residuals
r_espo <- espo_arma$residuals

par(mfrow=c(4,1))
plot(r_wti, main="WTI", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals")
plot(r_brent, main="BRENT", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals")
plot(r_dubai, main="Dubai", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals")
plot(r_espo, main="ESPO", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals")

plot(r_wti^2, main="WTI", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals^2")
plot(r_brent^2, main="BRENT", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals^2")
plot(r_dubai^2, main="Dubai", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals^2")
plot(r_espo^2, main="ESPO", 
     xlab="04/19/2013 - 03/06/2020",
     ylab="Residuals^2")

ArchTest(r_wti)
##WTI with ARCH effect

ArchTest(r_brent)
##Brent with no ARCH effect

ArchTest(r_dubai)
##Dubai with ARCH effect

ArchTest(r_espo)
##ESPO with ARCH effect


#Fit Univariate GARCH Models
# For WTI and BRENT
ugarch.spec1 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                          mean.model=list(armaOrder=c(2,1)),
                          distribution.model = "sstd")

wti_garch <- ugarchfit(spec = ugarch.spec1, data = return.train$WTI)
wti_garch

brent_garch <- ugarchfit(spec = ugarch.spec1, data = return.train$BRENT)
brent_garch

# For Dubai
ugarch.spec2 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(0,1)),
                           distribution.model = "sstd")
dubai_garch <- ugarchfit(spec = ugarch.spec2, data = return.train$DUBAI)
dubai_garch

# For ESPO
ugarch.spec3 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(1,1)),
                           distribution.model = "sstd")
espo_garch <- ugarchfit(spec = ugarch.spec3, data = return.train$ESPO)
espo_garch



dcc.garch.spec = dccspec(uspec = multispec(c(ugarch.spec1, 
                                             ugarch.spec1, 
                                             ugarch.spec2, 
                                             ugarch.spec3)),
                           dccOrder = c(1,1),
                           distribution = "mvt")

dcc.garch.spec
dcc.fit = dccfit(dcc.garch.spec, data = return.train)
plot(dcc.fit)


dcc.garch.spec2 = dccspec(uspec = multispec(c(ugarch.spec1, 
                                             ugarch.spec2)),
                         dccOrder = c(1,1),
                         distribution = "mvt")
dcc.fit = dccfit(dcc.garch.spec2, data = return.train[, c(1,3)])
plot(dcc.fit)


