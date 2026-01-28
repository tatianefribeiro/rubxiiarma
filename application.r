###############################################################################
#Paper: Reflected Unit Burr XII autoregressive moving average models
#Tittle: Application of the RUBXII-ARMA(p,q) to the actual data
#Author: Tatiane F. Ribeiro
#Last update: January 26, 2026
###############################################################################

# Clear the memory
rm(list = objects())

# Directory 
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/RUBXII-ARMA_EES_Springer_oct2025/rubxiiarma_app")

# Required packages
library(gdata) #allows to read excel files.  
library(e1071) #allows to calculate asymetry and kustosis.
library(tidyverse)
library(forecast)

# Useful functions and its sources
# ßARMA - available sources at
# https://github.com/vscher/barma           (1 and 2)
# http://www.ufsm.br/bayer/boot-barma.zip   (3)
source("barma.fit.r")    #1
source("barma.r")        #2
source("best.barma.r")   #3

# KARMA: available sources at
# https://github.com/fabiobayer/KARMA
source("kum-mu-phi.r")
source("karma.fit.r")
source("karma.r")
source("best.karma.r")
#source("best.karma_cristiane.r")

# RUBXIIARMA
source("rubxiiarma.fit.R")
source("best.rubxiiarma.r")
source("forec_measures.R")

# BSARMA
#source("bsarma.r")
#######################
# Data processing
#######################
# Data set
data = read_csv("Simples_Energia_Armazenda_Mês_data2004.csv") 
#View(data)
#data = readxl::read_xls("dataenergy.xls")
#data1<- data%>%filter(`UF Residência`==1)
#View(data)

#data <- data1
vr <- data$val_eaconsimp4

# Transform variable vr (stored energy rate) for unit interval
y1 <- vr/100
round(summary(y1),4)

# Monthly data
month = 4
year = 2004

# Tau value (If tau = 0.5, the conditional median is modeled)
tau = 0.5 

# Convert data in time series object
Y <- ts(y1, start = c(year,month), frequency = 12)

# Sample size (all observations)
n_first <- length(Y)

# Forecasting horizon
h1 <- 10

# Taking off the last h1 observations.
n <- n_first-h1    # Sample size for estimation

# Training data
y <- ts(y1[1:n],start = c(year,month),frequency = 12)  

# Validation data (to forecast)
yh <- ts(Y[(n+1):n_first],start = c(end(y)[1],end(y)[2]+1),
         frequency=12) 

######################################
# Descriptive statistic from data (y)
######################################
summary(y)  # resume
var(y)      # variance.
skewness(y) # skewness.
kurtosis(y) # kurtosis.

# Some graphics
hist(y)         
monthplot(y)     
# ggseasonplot(y, year.labels=TRUE, year.labels.left=TRUE) +
#   ylab("EAR(%)") +
#   ggtitle("Seasonal plot")

# TS plot
plot(y)

# Autocorrelation function.
Acf(y) 

# Partial autocorrelation function.
Pacf(y)

#######################################
#            Regressors
#######################################
# Time-varying regressors to model seasonality
t <- 1:length(y)       # in sample
t_hat <- (n+1):(n+h1)  # out of sample

C <- cos(2*pi*t/12)
C_hat <-cos(2*pi*t_hat/12)  

S <- sin(2*pi*t/12)
S_hat<-sin(2*pi*t_hat/12)  

# Regressors matrix
X <- cbind(C,S)
X_hat <- cbind(C_hat,S_hat)



#######################################
#           p and q orders
#######################################
# Choosing the best models from the RUBXII-ARMA, BARMA, and KARMA model classes
# pmax = 3
# qmax = 3
# rubxii_best <- best.rubxii(y, sf = c(start = c(year,month),frequency = 12),
#            h = h1, pmax = pmax, qmax = qmax,
#            nbest=10, tau=tau, link = "logit",
#            X = X, X_hat=X_hat)
# 
# barma_best <- best.barma(y, sf = c(start = c(year,month),frequency = 12),
#            h = h1, pmax = pmax, qmax = qmax,
#            nbest = 10, link = "logit",
#            X = X, X_hat=X_hat)
# 
# karma_best <- best.karma(y, sf = c(start = c(year,month),frequency = 12),
#            h = h1, pmax = pmax, qmax = qmax,
#            nbest = 10,link = "logit",
#            X = X, X_hat=X_hat)


# Orders of the final models 
p_rubxiiarma <- 1:3  #ordem 5: 3,5
q_rubxiiarma <- NA
p_karma <- 1:3
q_karma <- NA#1:3
p_barma <- 1:3#1:2
q_barma <- NA#1:3

# Final models
fit_rubxiiarma <- rubxiiarma.fit(y,ar = p_rubxiiarma,
                               ma = q_rubxiiarma,
                               tau = tau,
                               link = "logit", h = h1, 
                               diag = 2, #2:save the plots
                               X = X, X_hat=X_hat)

fit_barma <- barma(y,ar = p_barma, ma = q_barma, 
                   link = "logit",
                   h = h1, diag = 2,
                   X = X, X_hat = X_hat)

fit_karma <- karma(y, ar = p_karma, ma = q_karma, link = "logit", 
                   h = h1, diag = 2, 
                  X = X, X_hat = X_hat)



# Diagnostic measures (quantile residuals)
res_RUBXII = fit_rubxiiarma$residuals
(bltest_RUBXII <- Box.test(res_RUBXII, lag = 20, type =  "Ljung-Box"))#, fitdf = 2)

res_KARMA = fit_karma$resid3 
(bltest_KARMA <- Box.test(res_KARMA, lag = 20, type =  "Ljung-Box"))#, fitdf = 2)

res_BARMA = fit_barma$resid5  
(bltest_BARMA <- Box.test(res_BARMA, lag = 20, type =  "Ljung-Box"))#, fitdf = 2)


Acf(res_RUBXII)
Pacf(res_RUBXII)
Acf(res_BARMA)
Pacf(res_BARMA)
Acf(res_KARMA)
Pacf(res_KARMA)

fit_rubxiiarma$aic
fit_karma$aic
fit_barma$aic

fit_rubxiiarma$bic
fit_karma$bic
fit_barma$bic
################
# Tables Latex
###############
# Fits
round(rbind(
fit_rubxiiarma$model,
fit_karma$model,
fit_barma$model),4)

stargazer::stargazer(
  rbind(
  fit_rubxiiarma$model,
  fit_karma$model,
  fit_barma$model),digits=4)

stargazer::stargazer(fit_bsarma$model, digits=4)
# AIC, BIC and Box-Ljung test p-value 
diag_meas <- matrix(NA,3,3)
diag_meas<- cbind(
  rbind(
    fit_rubxiiarma$aic,
    fit_karma$aic,
    fit_barma$aic),
  rbind(fit_rubxiiarma$bic,
        fit_karma$bic,
        fit_barma$bic),
  rbind(bltest_RUBXII$p.value,
        bltest_KARMA$p.value,
        bltest_BARMA$p.value)) 
colnames(diag_meas) <- c("AIC","BIC","pvalue")
rownames(diag_meas) <- c("RUBXII-ARMA", "KARMA", "BARMA")
stargazer::stargazer(diag_meas,digits=4)

# Forecasting adequacy measures
acc_fits <- round(
  rbind(accuracy(fit_rubxiiarma$forecast,yh)[,c(2:3,5)],
        accuracy(fit_karma$forecast,yh)[,c(2:3,5)],
        accuracy(fit_barma$forecast,yh)[,c(2:3,5)]),4)

acc_fits
stargazer::stargazer(acc_fits,digits=4)

et=yh-fit_rubxiiarma$forecast
pt=100*et/yh

(rmse=sqrt(1/h1*sum(et^2)))
(mae=1/h1*sum(abs(et)))
(mpe=1/h1*sum(pt))
(mape=1/h1*sum(abs(pt)))


##############
# PLOTS
#############
# TIME SERIES PLOT 
w1<-9
h2<-7
postscript(file = "y_plots.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(2,2))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  plot(y,ylab = "Stored hydroelectric energy")
  monthplot(y,ylab = "Stored hydroelectric energy")
  Acf(y,lag=30, main="")
  Pacf(y, main="")
}
dev.off()

# Sampling ACF and PACF of the residuals
w1<-11
h2<-5
postscript(file = "residuals_plots.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,2))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  Acf(res_RUBXII)
  Pacf(res_RUBXII)
  # Acf(res_KARMA)
  # Pacf(res_KARMA)
  # Acf(res_BARMA)
  # Pacf(res_BARMA)
}
dev.off()


w1<-8
h2<-5
postscript(file = "fitt_value_foresc.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  plot(y,type="l",ylab="Stored hydroelectric energy",xlab="Time",
       ylim=c(min(y),max(y)))
  lines(fit_rubxiiarma$fitted,col="blue",lty=2)
  fim<-end(y)[1]+end(y)[2]/12
  abline(v=fim,lty=2)
  lines(yh,col="blue",lty=2)
  legend("bottomleft",c("Observed data","Fitted values"),#pch=vpch,
         pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue",
                                                  "red"))
}
dev.off()


w1<-7
h2<-6
postscript(file = "resid_vs_index.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  t<-seq(-5,n+h1-6,by=1)
  plot(res_RUBXII,main=" ",xlab="Index",ylab="Quantile residuals",
       pch = "+",ylim=c(-5,5))
  lines(t,rep(-3,n+h1),lty=2,col=1)
  lines(t,rep(3,n+h1),lty=2,col=1)
  lines(t,rep(-2,n+h1),lty=3,col=1)
  lines(t,rep(2,n+h1),lty=3,col=1)

}
dev.off()


w1<-11
h2<-8
postscript(file = "fitted_values_res_plots.eps",horizontal=F,
           paper="special",width = w1, height = h2,
           family = "Times",
           pointsize = 15)
{
  par(mfrow=c(2,2))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  #
  plot(y,type="l",ylab="Stored hydroelectric energy",xlab="Time",
       ylim=c(min(y),max(y)))
  lines(fit_rubxiiarma$fitted,col="blue",lty=2)
  fim<-end(y)[1]+end(y)[2]/12
  abline(v=fim,lty=2)
  lines(yh,col="blue",lty=2)
  legend("bottomleft",c("Observed data","Fitted values"),#pch=vpch,
         pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue",
                                                  "red"))
  #
  t<-seq(-5,n+h1-6,by=1)
  plot(res_RUBXII,main=" ",xlab="Index",ylab="Quantile residuals",
       pch = "+",ylim=c(-5,5))
  lines(t,rep(-3,n+h1),lty=2,col=1)
  lines(t,rep(3,n+h1),lty=2,col=1)
  lines(t,rep(-2,n+h1),lty=3,col=1)
  lines(t,rep(2,n+h1),lty=3,col=1)
  #
  Acf(res_RUBXII)
  Pacf(res_RUBXII)
}
dev.off()

