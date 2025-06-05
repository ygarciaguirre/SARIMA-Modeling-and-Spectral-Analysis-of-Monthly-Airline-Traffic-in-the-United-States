library(astsa)
library(MASS)
library(forecast)

# Accessing the data set
View(air_traffic)
air_traf <- air_traffic$Pax[73:192]/1000000         # Values in millions
air <- ts(air_traf, start=c(2009,1), frequency=12)  # Time series from 2009-2018
air

# Box-Jenkins approach
# 1. Plots
plot.ts(air)
acf(air)

decomposition <- decompose(air)
plot(decomposition)

decomposition
plot(decomposition$x)
plot(decomposition$seasonal)
plot(decomposition$trend)

# 2. Transformations
transform <- boxcox(air~as.numeric(1:length(air)))
opt_lambda <- transform$x[which(transform$y == max(transform$y))]
opt_lambda                             # optimal lambda = 0.87

bc_air <- (air^opt_lambda-1)/opt_lambda
plot.ts(bc_air)                        # the transformation decreases volatility

# library(forecast)
# library(zoo)
trend <- rollmean(bc_air, k=12, fill=NA, align='right')
air_trend <- autoplot(trend)
plot(air_trend)                        # show an upward trend

air_season <- ggseasonplot(bc_air)
plot(air_season)                       # shows a monthly seasonality

air_d1 <- diff(bc_air, 1) 
plot.ts(air_d1)                        # gets rid of the trend

air_d12 <- diff(air_d1, 12)
plot.ts(air_d12)                       # makes it even better

Box.test(air_d12, type="Ljung-Box") # test stationary signal

var(air); var(bc_air); var(air_d1); var(air_d12) # we have decreased the variance

plot.ts(air_d12)
hist(as.numeric(air_d12), breaks=seq(from=-3, to=3, by=1/5))

variances <- c(round(var(air),4), var(bc_air), var(air_d1), var(air_d12))
variances <- round(variances, 4)
variances

acf(bc_air); acf(air_d1); acf(air_d12)

pacf(bc_air); pacf(air_d1); pacf(air_d12)


# 3. Identification of orders of the model
# We use a SARIMA model since we differenced once to remove trend
# and then a second time to remove seasonality --> d=1 and D=1

acf2(air_d12)

# Seasonal component
  # ACF cuts off after h=1 --> SMA component with Q=1
  # PACF has no significant spikes --> no SRA component --> P=0

# Non-seasonal component
  # ACF cuts off after lag s=1 --> MA component with q=1
  # PACF is tailing off at s=1,2,3,... --> no AR component --> p=0
      # OR
  # ACF is tailing off at s=1,4,5,6,... --> no MA component --> q=0
  # PACF cuts off after lag s=3 --> AR component with p=3
      # OR
  # Both tailing off --> p=1 and q=1

auto.arima(bc_air)
  # ARIMA(0,1,1)(1,1,2)[12]

# Possible models
  # Model 1: ARIMA(0,1,1)(0,1,1)[12]
  # Model 2: ARIMA(3,1,0)(0,1,1)[12]
  # Model 3: ARIMA(0,1,1)(1,1,2)[12]
  # Model 4: ARIMA(1,1,1)(0,1,1)[12]
  # Model 5: ARIMA(0,1,2)(0,1,1)[12]
  # Model 6: ARIMA(3,1,1)(0,1,1)[12]
  # Model 7: ARIMA(4,1,0)(0,1,1)[12]
  # Model 8: ARIMA(2,1,0)(0,1,1)[12]

# 4. Estimation of parameters
# ARIMA(0,1,1)(0,1,1)[12]
model_1 <- sarima(bc_air, p=0, d=1, q=1, P=0, D=1, Q=1, S=12)
model_1
# both coefficients are significant
# AIC = 1.844837  AICc = 1.845915  BIC = 1.919776

# ARIMA(3,1,0)(0,1,1)[12]
model_2 <- sarima(bc_air, p=3, d=1, q=0, P=0, D=1, Q=1, S=12)
model_2
# all coefficients are significant
# AIC = 1.86264  AICc = 1.866305  BIC = 1.987538

# ARIMA(0,1,1)(1,1,2)[12]
model_3 <- sarima(bc_air, p=0, d=1, q=1, P=1, D=1, Q=2, S=12)
model_3
# SAR1 and SMA1 are not significant
# AIC = 1.822127  AICc = 1.825792  BIC = 1.947025

# ARIMA(1,1,1)(0,1,1)[12]
model_4 <- sarima(bc_air, p=1, d=1, q=1, P=0, D=1, Q=1, S=12)
model_4
# AR1 is not significant
# AIC = 1.860452  AICc = 1.86263  BIC = 1.960371

# ARIMA(0,1,2)(0,1,1)[12]
model_5 <- sarima(bc_air, p=0, d=1, q=2, P=0, D=1, Q=1, S=12)
model_5
# MA2 is not significant
# AIC = 1.859544  AICc = 1.861721  BIC = 1.959462

# ARIMA(3,1,1)(0,1,1)[12]
model_6 <- sarima(bc_air, p=3, d=1, q=1, P=0, D=1, Q=1, S=12)
model_6
# AR3 and MA1 are not significant
# AIC = 1.881031  AICc = 1.886583  BIC = 2.030909

# ARIMA(4,1,0)(0,1,1)[12]
model_7 <- sarima(bc_air, p=4, d=1, q=0, P=0, D=1, Q=1, S=12)
model_7
# AR3 and AR4 are not significant
# AIC = 1.881073  AICc = 1.886625  BIC = 2.030951

# ARIMA(2,1,0)(0,1,1)[12]
model_8 <- sarima(bc_air, p=2, d=1, q=0, P=0, D=1, Q=1, S=12)
model_8
# all components are significant
# AIC = 1.903799  AICc = 1.905977  BIC = 2.003718

model_1$ICs; model_2$ICs; model_3$ICs; model_4$ICs; model_5$ICs; model_6$ICs; model_7$ICs; model_8$ICs

# AIC and AICc both prefer model 4, but it does not pass the coefficient test
# BIC prefers model 1
# We continue with models 1 and 2

# 5. Diagnostics
sarima(bc_air, p=0, d=1, q=1, P=0, D=1, Q=1, S=12)
# Plot
  # There are no obvious patterns, except for the data between 2008 and 2010
  # There are two outliers by the end of 2013 and the end of 2017 that exceed 3 standard deviations
  # However, there are no values that exceed in magnitude
# ACF
  # It agrees with our model assumptions
# Q-Q plot
  # Since most data points lie on the normal lie, we can assume normality despite the two outliers
# p-values for the Q test
  # p-values for lags between 3 and 36
  # They all exceed 0.05, so we can assume that the residuals are white noise

sarima(bc_air, p=3, d=1, q=0, P=0, D=1, Q=1, S=12)
# Similar plots, but we have
  # ACF has two borderline spikes
  # There is a borderline p-value, and most are smaller

# 6. Model selection
model_1$ICs; model_2$ICs

# Comparing the AIC, AICc, and BIC values, they all prefer model 1
# Also, the residual diagnostics prefers model 1
# We conclude that our data set follows a ARIMA(0,1,1)(0,1,1)[12] process

# Forecast
# ARIMA(0,1,1)(0,1,1)[12]
sarima.for(bc_air, n.ahead=12, p=0, d=1, q=1, P=0, D=1, Q=1, S=12)
sarima.for(bc_air, n.ahead=12, p=0, d=1, q=1, P=0, D=1, Q=1, S=12, plot.all=TRUE, gg=FALSE)

forecast <- sarima.for(bc_air, n.ahead=12, p=0, d=1, q=1, P=0, D=1, Q=1, S=12, plot=FALSE)
ts.plot(bc_air, forecast$pred, col=1:2, xlim=c(2009,2020))
U1 <- forecast$pred + 2*forecast$se  # 95% confidence interval
L1 <- forecast$pred - 2*forecast$se  # 2 sd ~ 1.96
xx1 <- c(time(U1), rev(time(U1)))
yy1 <- c(L1, rev(U1))
polygon(xx1, yy1, border = 8, col = gray(.6, alpha = .2))
lines(forecast$pred, type="p", pch=20, col=2)

air_test <- air_traffic$Pax[193:204]/1000000        # 2019 data for forecasting
air_test <- ts(air_test, start=c(2019,1), frequency=12)
bc_air_test <- (air_test^opt_lambda-1)/opt_lambda

ts.plot(bc_air, forecast$pred, col=1:2, xlim=c(2018,2020))
polygon(xx1, yy1, border = 8, col = gray(.6, alpha = .2))
lines(forecast$pred, type="p", pch=20, col=2)
points(bc_air_test, pch=20, col="blue")


# Spectral Analysis
acf(bc_air)

library(TSA)
my_periodogram <- periodogram(bc_air, log='no', plot=TRUE, ylab="Periodogram",
                              xlab="Frequency", lwd=3)
my_periodogram
first_peak <- which.max(my_periodogram$spec)
first_peak
second_peak <- which.max(my_periodogram$spec[first_peak])
second_peak

1/my_periodogram$freq[first_peak]
1/my_periodogram$freq[second_peak]

0.1/12

1/(2*(0.1/12))
1/(11*(0.1/12))

air_per <- mvspec(bc_air, taper=0.1, log="no")
spectrum(bc_air, log="no")
air_smo <- mvspec(bc_air, taper=0.1, spans=c(5,5), log="no")
spectrum(bc_air, spans=c(5,5), log="no")

# cex.main=0.5

log_air_per <- mvspec(bc_air, taper=0.1, log="yes")
spectrum(bc_air)
log_air_smo <- mvspec(bc_air, taper=0.1, spans=c(5,5), log="yes")
spectrum(bc_air, spans=c(5,5))
