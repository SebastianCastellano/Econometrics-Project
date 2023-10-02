# Code for the Final Project of Econometrics course 2021/2022
# Group 8: Sebastian Castellano, Giulia Mulattieri, Virginia Muscionico

#Effects of the US Stock Market Return and Volatility on the VKOSPI

# 0) IMPORT LIBRARIES
rm(list=ls())
graphics.off()
library(tidyverse)
library(pracma)
library(zoo)
library(readxl)
library(lmtest)
library(tseries)
library(urca) 
library(moments)
library("PerformanceAnalytics")
library(corrplot)
library(strucchange)
library(dynlm)
library(orcutt)
library(sandwich)
library(xts)
library(quantmod)
library(fUnitRoots)
library(FinTS)
library(rugarch)
library(ggplot2)
library(forecast)
library(fGarch)
library(alfred)


#1) IMPORT DATA and create dataset ####
data <- read_excel('Vkospi_data.xlsx')
InflationUSAKOR.monthly <- read_csv("InflationUSAKOR.csv")
InflationUSA.monthly <- InflationUSAKOR.monthly %>% dplyr::filter(LOCATION == 'USA')
InflationKOR.monthly <- InflationUSAKOR.monthly %>% dplyr::filter(LOCATION == 'KOR')

data %>% left_join(select(InflationUSA.monthly,TIME,Value),by="TIME") -> data

data %>% left_join(select(InflationKOR.monthly,TIME,Value),by="TIME") -> data

data %>%
  dplyr::rename(InflationUSA = Value.x) %>%
  dplyr::rename(InflationKOR = Value.y) -> data

data$Date <-as.Date(as.character(data[["Date"]]), "%Y%m%d")
Date <- data$Date
my_data_used <- data[c(-1,-6,-7,-8,-9,-14)]
# you should now have a data.frame called data with dimensions 2430*16

#2) OBSERVE DATA ####
#Plot of Vkospi and Vix
x11()
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(data$Date,data$VKOSPI, type='l', col="blue",xlab = "Date", ylab= "VKOSPI", main="VKOSPI and VIX")
par(new = TRUE)
plot(data$VIX, type = 'l', col="red", axes=FALSE, xlab = "", ylab = "")
axis(side=4)
mtext("VIX", side=4, line=3)

#Plot of SP and kospi
x11()
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(data$Date,data$SP500, type='l', col="black",xlab = "Date", ylab= "S&P500", main="S&P500 and KOSPI200")
par(new = TRUE)
plot(data$KOSPI200, type = 'l', col="darkgreen", axes=FALSE, xlab = "", ylab = "")
axis(side=4)
mtext("KOSPI200", side=4, line=3)

#Plot of sp500 and vix
x11()
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(data$Date,data$SP500, type='l', col="black",xlab = "Date", ylab= "S&P500", main="S&P500 and VIX")
par(new = TRUE)
plot(data$VIX, type = 'l', col="red", axes=FALSE, xlab = "", ylab = "")
axis(side=4)
mtext("VIX", side=4, line=3)

#Plot of Kospi200 and vkospi
x11()
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(data$Date,data$KOSPI200, type='l', col="black",xlab = "Date", ylab= "KOSPI200", main="KOSPI200 and VKOSPI")
par(new = TRUE)
plot(data$VKOSPI, type = 'l', col="blue", axes=FALSE, xlab = "", ylab = "")
axis(side=4)
mtext("VKOSPI", side=4, line=3)

# Correlation plots
x11()
corrplot(cor(my_data_used))
x11()
chart.Correlation(my_data_used, histogram=TRUE, pch=19, col="blue")


# 3) DATA CLEANING and Hypothesis checking ####
#a) there is a weird outlier in STL, probably human error since it is 10 times all others. We remove it
x11()
boxplot(data$STL, col="gold")

x11()
par(mfcol=c(1,2))
plot(Date,data$STL, col=ifelse(data$STL>100,"red","black"), pch=18, main="STL before")
stl_nooutliers<- data$STL[252] / 10
data$STL[252] <- stl_nooutliers
plot(Date,data$STL, pch=18, main="STL after")


#b) Check logreturnExchangerate
#These two values below are very high with respect to all others but they are not errors.
#By checking the historic exchange rate we find the same value
x11()
plot(data$logreturnExchangerate)
data$logreturnExchangerate[1131]
data$logreturnExchangerate[1141]


#c) Check for structural changes (non stationarity) - Chow test for constant mean vs changing mean through time

# Plot the F-stat for the Chow test on vkospi
y <- data$VKOSPI
x11()
par(mfcol=c(2,1))
plot(Date,y, main="VKOSPI", col="blue",type="l")
y_ts<-ts(y)
y.ftest <- Fstats(y_ts ~ 1)
plot(y.ftest, main="F statistic for Chow test") #when the graph is above the red line, we reject the hypothesis of stationarity
#plot(y.ftest, pval = TRUE, ylim = c(0,0.1)) #complementary to above, we reject if below

adf.test(y_ts)$p.value

# Check for all the timeseries
for (i in 1:dim(my_data_used)[2]){
  y <- my_data_used[[i]]
  print(names(my_data_used)[[i]])
  y.ftest <- Fstats(y ~ 1)
  print(sctest(y.ftest, type = "aveF"))
  print(adf.test(y)$p.value)
}
#we observe that we reject stationarity for all, with the exception of logreturn(Exchangerate)
# Makes sense since it is already transformed


#4) REGRESSORS AND VARIABLES ####
#a) Create our regressors following the paper
graphics.off()
y <- log(data$VKOSPI)
rek <- c(NA,diff(log(data$KOSPI200), lag = 1))
reus <- c(NA,diff(log(data$SP500), lag = 1))
lvix <- log(data$VIX)
ex <- data$logreturnExchangerate
cdrate <- data$`CD(91D)`
termspread <- data$`termspread(5Y-91D)`
creditspread <- data$`creditspread(BBB-AA)`
opi<-data$OPI
stl<-data$STL
volume<-data$Volume
yZ <- data.frame(y=y,
                 rek=rek,
                 reus=reus,
                 lvix=lvix,
                 ex=ex,
                 cdrate=cdrate,
                 termspread=termspread
                 )

corrplot(cor(na.omit(yZ)))

#b) Chow Test and ADF test on the transformed timeseries
x11()
par(mfcol=c(2,1))
plot(Date,y, main="VKOSPI", col="blue",type="l")
y_ts<-ts(y)
y.ftest <- Fstats(y_ts ~ 1)
plot(y.ftest, main="F statistic for Chow test") 
adf.test(y_ts)$p.value

for (i in 1:dim(yZ)[2]){
  y <- na.omit(yZ[[i]])
  print(names(yZ)[[i]])
  y.ftest <- Fstats(y ~ 1)
  print(sctest(y.ftest, type = "aveF"))
  print(adf.test(y)$p.value)
}
# Now the diff-log transformed variables are stationary at least

#c) Compute simple moving averages
y5 <- rollmean(y,5, fill = NA, align = 'right')
y10 <- rollmean(y,10, fill = NA, align = 'right')
y22 <- rollmean(y,22, fill = NA, align = 'right')

#d) Create other macro-economic variables
inflationKOR.USA <- data$InflationKOR/data$InflationUSA


#5) MODELS ####
#a) Models from the paper
X1 <- cbind(y, y10)  #1st model from matlab files
X2 <- cbind(y, y10, ex, cdrate, termspread, creditspread)  #2nd model from matlab files
X3 <- cbind(y, lvix, y5, ex)  #3rd model from matlab files
X4 <- cbind(y, reus, y10, cdrate, creditspread)  #4th model from matlab files
X5 <- cbind(y, rek, y10, ex, cdrate, termspread, creditspread)  #5th model from matlab files
X6 <- cbind(y, reus, lvix, y10, cdrate)  #6th model from matlab files
X7 <- cbind(y, rek, lvix, y5, ex)  #7th model from matlab files

#b) Our additional models
Our1 <- cbind(rek, reus, lvix, ex, cdrate, termspread, creditspread)
Our2 <- cbind(rek, reus, lvix, ex, cdrate, creditspread)
Our3 <- cbind(rek,reus,lvix,ex,cdrate,termspread,creditspread,inflationKOR.USA,stl,volume)
Our4 <-cbind(rek,reus,lvix,ex,termspread,creditspread,inflationKOR.USA,stl)
Our5 <- cbind(rek,reus,lvix,cdrate, inflationKOR.USA)
Our6 <-cbind(lvix,termspread,inflationKOR.USA,creditspread,reus)

#c) Regression
#Select hyperparameters
time_window = 1008
steps_ahead = 5
X <- X6  #<----------  CHANGE THIS TO SELECT THE MODEL YOU WANT TO USE

#Analysis of lm for 1-step ahead, full training set
mod <- lm(y[2:length(y)] ~ X[-dim(X)[1],], na.action = "na.exclude")
summary(mod)

x11()
par(mfrow=c(2,2))
plot(mod)

#6) TESTS ####
# Breusch-Pagan Test:
bptest(mod)$p.value #(H_0: homoschedasticity)

# Box test
Box.test(residuals(mod))$p.value #(H_0: independence of the residuals)

#Shapiro test
shapiro.test(residuals(mod))$p.value #(H_0: normality)


#7) PREDICT ####
# Create X and Y matrices
X <- X[-((dim(X)[1]-steps_ahead +1):dim(X)[1]),]  #remove the last k=steps_ahead values
Y <- y[-(1:steps_ahead)]  #remove the first k=steps_ahead values
reg <- lm(Y ~ X, na.action = "na.exclude")  #na.exclude gives same number of row in output
summary(reg)

# Predict
Y.pred = rep(0,(length(Y)-time_window))
for (i in 1:(length(Y)-time_window)){
  X_window <- X[(i:(i+time_window-1)),]
  Y_window <- Y[(i:(i+time_window-1))]
  reg <- lm(Y_window ~ X_window, na.action = "na.exclude")
  Y.pred[i] <- as.numeric(X[i+time_window,]%*%reg$coefficients[-1]+reg$coefficients[1])
}

Y.actual <- Y[(time_window+1):length(Y)]
mse <- mean((Y.actual-Y.pred)^2)
mae <- mean(abs(Y.actual-Y.pred))
print(mae)
print(mse)

x11()
plot(Y.actual, type = 'l',col= "green",main = paste("Actual vs Predicted - Steps Ahead =",steps_ahead))
lines(Y.pred, col="red")


#8) GARCH and ARCH ####
# Chow test - Fstats
d_y_ts=diff(y_ts) #diff log of the ratio
x11()
plot(d_y_ts,ylab="logret",col="gold",lwd=2, main="Diff(Log(VKOSPI))")
abline(h=0,lty=2,col="red")

x11()
fs.d_y<-Fstats(d_y_ts~1) 
plot(fs.d_y) #now the plot is always under the red line -> stationarity!

# Structural change 
sctest(fs.y)#low pvalue reject the null hypotesis -> there is structural change in our time series
plot(y,type="l",col="blue")
lines(breakpoints(fs.y)) #takes the highest value of fs.y -> break point
breakpoints(fs.y)

#Diff(Log(VKOSPI)) with breakline
plot(d_y_ts,col="gold", main="Diff(Log(VKOSPI))")
lines(breakpoints(fs.d_y))
abline(h=0,lty=2,col="red")
sctest(fs.d_y,type="aveF")

#ACF
# We first do for log(vkospi) out of curiosity but remember that we do not have stationarity
x11()
par(mfrow=c(1,2))
plot(acf(y_ts,12,xlim=c(1,13)),main="ACF of log(y)")
plot(acf(y_ts,12,type="partial",xlim=c(1,13))) #first time step is partially autocorrelated with itself

res<-acf(y_ts,12,type="partial",plot=FALSE)
res

# For dif(log(vkospi)) we have stationarity
x11()
par(mfrow=c(1,2))
plot(acf(d_y_ts,12,xlim=c(1,13))) # 1 at 0 and a significant one at time step 5  
plot(acf(d_y_ts,12,type="partial",xlim=c(1,13)))

adf.test(d_y_ts) #stationary time series: low pvalue
adf.test(y_ts) #non stationary: high pvalue

# Check for arch
archtest<-ArchTest(d_y_ts,lags=1,demean=TRUE)
archtest #there is an arch effect in d_y_ts

#a) Estimate mean equation r=beta+error
mod<-dynlm(d_y_ts~1)
#b) Retrieve the residuals from the former model
ehatsq<-ts(resid(mod)^2)
#c) Regress squared residuals on one-lagged squared
mod_arch<-dynlm(ehatsq~L(ehatsq),data=ehatsq)
summary(mod_arch)

# Stat significant 
# R2 squared very low and high volatility
x11()
arch.fit=garchFit(~arma(0,0)+garch(1,1),data=d_y_ts,trace=F)
summary(arch.fit)
pred.arch <- predict(arch.fit,n.ahead=20,plot=TRUE)
pred.mod_arch <- predict(mod_arch,n.ahead=20,plot=TRUE)
