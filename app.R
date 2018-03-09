# Antonis Pishias
# Financial Data and Statistics
# Assignment 2

#
# Section # Importing ##################################
#

setwd("~/Sites/rstudio/fds2")
setwd("./data//cryptocurrencyCSV/")

# install.packages("xts")
# install.packages("devtools")
# library(devtools)
# install_github("cran/PerformanceAnalytics")
# install.packages("astsa")
# install.packages("pracma")
# install.packages("corrplot")
# install.packages("fractaldim")
# install.packages("forecast")
# install.packages("fpp2")
library(xts)
library(astsa)
# library(PerformanceAnalytics)
# library(pracma)
library(corrplot)
library(fractaldim)
library(forecast)
library(fpp2)

files <- list.files()
cryptos_daily_lr <- list()
cryptos_daily_close <- list()
cryptos_weekly_lr <- list()
cryptos_weekly_close <- list()
cryptos_monthly_lr <- list()
cryptos_monthly_close <- list()

# sample <- read.csv(files[1], nrows=5)
# as.Date(sample$Date, format = "%b %d ,%Y")

for (file in files){
  crypto <- gsub("*.csv", "", file)
  zooObj <- read.zoo(file, index.column = "Date", sep = ",", header = T,
                     format="%b %d, %Y")
  zooObj <- na.omit(zooObj$Close)
  daily_close <- as.xts(as.numeric(zooObj), order.by = index(zooObj))
  colnames(daily_close) <- "Close"
  
  # Daily Log-Returns
  cryptos_daily_close[[crypto]] <- daily_close
  xts_dlr <- diff(log(daily_close$Close))
  colnames(xts_dlr) <- c("Daily Log-Returns")
  head(xts_dlr, n=5)
  cryptos_daily_lr[[crypto]] <- xts_dlr
  
  # Weekly Log-Returns
  weekly_close <- to.weekly(daily_close, OHLC = FALSE)
  cryptos_weekly_close[[crypto]] <- weekly_close
  
  head(weekly_close, n=5)
  
  xts_wlr <- diff(log(weekly_close$Close))
  colnames(xts_wlr) <- c("Weekly Log-Returns")
  cryptos_weekly_lr[[crypto]] <- xts_wlr
  
  # Monthly Log-Returns 
  monthly_close <- to.monthly(daily_close, OHLC = FALSE)
  cryptos_monthly_close[[crypto]] <- monthly_close
  
  head(monthly_close, n=5)
  
  xts_mlr <- diff(log(monthly_close$Close))
  colnames(xts_mlr) <- c("Monthly Log-Returns")
  cryptos_monthly_lr[[crypto]] <- xts_mlr
  
  rm(daily_close, weekly_close, monthly_close, xts_dlr, xts_wlr, xts_mlr, zooObj, crypto)
}
rm(file, files)
setwd("./../")

#
# Section # Helper Functions ###########################
#

mergeAll <- function(crypto_xts_list) {
  m1 <- merge(crypto_xts_list$Bitcoin, crypto_xts_list$Dash, join = "left")
  m1 <- merge(m1, crypto_xts_list$Ethereum, join = "inner")
  m1 <- merge(m1, crypto_xts_list$Litecoin, join = "inner")
  m1 <- merge(m1, crypto_xts_list$Monero, join = "inner")
  m1 <- merge(m1, crypto_xts_list$NEM, join = "inner")
  m1 <- merge(m1, crypto_xts_list$Ripple, join = "inner")
  m1 <- merge(m1, crypto_xts_list$Siacon, join = "inner")
  m1 <- merge(m1, crypto_xts_list$Stellar, join = "inner")
  m1 <- merge(m1, crypto_xts_list$Verge, join = "inner")
  colnames(m1) <- names(crypto_xts_list)
  m1 <- na.omit(m1)
  m1
}

## put histograms on the diagonal
panel.hist <- function(x)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (r+0.7))
}

#
# Section # Time-Horizons ##############################
#

# Daily #
daily_lr <- mergeAll(cryptos_daily_lr)
periodicity(daily_lr)

# Scatter Plots
dfc <- as.data.frame(daily_lr)
pairs(dfc,diag.panel = panel.hist, upper.panel = panel.cor)

# Time-series Plots
dfcTS <- ts(dfc, frequency = 365, start = c(2015,06,26))
# Pairs
autoplot(dfcTS, facets=T,main= "Cryptocurrency Daily Log-Returns") + ylab("Log-Returns") + xlab("Time")

# Weekly #
weekly_lr <- mergeAll(cryptos_weekly_lr)
periodicity(weekly_lr)

# Scatter Plots
wfc <- as.data.frame(weekly_lr)
pairs(wfc,diag.panel = panel.hist, upper.panel = panel.cor)

# Time-series Plots
wfcTS <- ts(wfc, frequency = 52, start = c(2015,06,26))
# Pairs
autoplot(wfcTS, facets=T,main= "Cryptocurrency Weekly Log-Returns") + ylab("Log-Returns") + xlab("Time")

# Monthly #
monthly_lr <- mergeAll(cryptos_monthly_lr)
periodicity(monthly_lr)

# Scatter Plots
mfc <- as.data.frame(monthly_lr)
pairs(mfc,diag.panel = panel.hist, upper.panel = panel.cor)

# Time-series Plots
mfcTS <- ts(mfc, frequency = 12, start = c(2015,06,26))
# Pairs
autoplot(mfcTS, facets=T,main= "Cryptocurrency Monthly Log-Returns")  + ylab("Log-Returns") + xlab("Time")

#
# Section # Export #####################################
#

daily_cl <- as.data.frame((mergeAll(cryptos_daily_close)))
weekly_cl <- as.data.frame((mergeAll(cryptos_weekly_close)))
monthly_cl <- as.data.frame((mergeAll(cryptos_monthly_close)))

setwd("./output/")
write.csv(daily_cl, "cryptos_daily_close.csv", row.names = TRUE)
write.csv(weekly_cl, "cryptos_weekly_close.csv", row.names = TRUE)
write.csv(monthly_cl, "cryptos_monthly_close.csv", row.names = TRUE)

write.csv(dfc, "cryptos_daily_lr.csv", row.names = TRUE)
write.csv(wfc, "cryptos_weekly_lr.csv", row.names = TRUE)
write.csv(mfc, "cryptos_monthly_lr.csv", row.names = TRUE)
setwd("./../")

#
# Section # Auto-Correlerations ########################
#

NEM_log_returns <- daily_lr$NEM["2015-06-26/"]
acf2(NEM_log_returns, max.lag = 60)

plot(abs(NEM_log_returns))
acf2(abs(NEM_log_returns), max.lag = 60)

plot(NEM_log_returns^2)
acf2(NEM_log_returns^2, max.lag = 60)

# Ljung-Box test 
# a p-value greater than 0.05 suggests that the data are not significantly different from white noise.
Box.test(NEM_log_returns, lag = 10, fitdf = 0, type = "Ljung")

#
# Section # ARIMA ######################################
#

# This section assumes normality
# See 'arimaplus.m' for ARIMA with t-student residuals

# NEM Daily
test = 30
d_cl <- ts(log(cryptos_daily_close$NEM["2015-06-26/"]), frequency = 365, start = c(2015,06,26))
size <- length(d_cl)
d_cl_train <- ts(d_cl[1:(size-test)], frequency = 365, start = c(2015,06,26))
str(d_cl_train)
acf2(d_cl_train, max.lag = 60)

d1_model <- sarima(d_cl_train, 1,1,0)
d1_model$ttable
d1_model$AIC
d1_model$BIC

sarima.for(d_cl_train, n.ahead=test, 1,1,0)
lines(d_cl)

# NEM Weekly
test = 15
w_cl <- ts(log(cryptos_weekly_close$NEM["2015-06-26/"]), frequency = 52, start = c(2015,06,26))
size <- length(w_cl)
w_cl_train <- ts(w_cl[1:(size-test)], frequency = 52, start = c(2015,06,26))
str(w_cl_train)
acf2(w_cl_train, max.lag = 60)

w1_model <- sarima(w_cl_train, 1,1,0)
w1_model$ttable
w1_model$AIC
w1_model$BIC

sarima.for(w_cl_train, n.ahead=test, 1,1,0)
lines(w_cl)

# NEM Monthly
test = 15
m_cl <- ts(log(cryptos_monthly_close$NEM["2015-06-26/"]), frequency = 12, start = c(2015,06,26))
size <- length(m_cl)
m_cl_train <- ts(m_cl[1:(size-test)], frequency = 12, start = c(2015,06,26))
str(m_cl_train)
acf2(m_cl_train)

m1_model <- sarima(m_cl_train, 1,1,0)
m1_model$ttable
m1_model$AIC
m1_model$BIC

sarima.for(m_cl_train, n.ahead=test, 1,1,0)
lines(m_cl)

#
# Section # Correlations ###############################
#

# Daily #
dfc_cor <- cor(dfc, use = "everything", method = "pearson")
corrplot(dfc_cor, method="color", type = "upper", title="Daily Log-Return Correlation")

# Weekly #
wfc_cor <- cor(wfc, use = "everything", method = "pearson")
corrplot(wfc_cor, method="color", type = "upper", title="Weekly Log-Return Correlation")

# Monthly #
mfc_cor <- cor(mfc, use = "everything", method = "pearson")
corrplot(mfc_cor, method="color", type = "upper", title="Monthly Log-Return Correlation")

#
# Section # Hurst ######################################
#

hurstexp(daily_lr$NEM)
hurstexp(weekly_lr$NEM)
hurstexp(monthly_lr$NEM)

#
# Section # Scaling # Log-returns relative frequency ###
#

d_histinfo <- hist(daily_lr$NEM, breaks=170)
w_histinfo <- hist(weekly_lr$NEM, breaks=25)
m_histinfo <- hist(monthly_lr$NEM, breaks=5)
l_rnd <- rlaplace(100000, location=0, scale=.05)
l_histinfo <- hist(l_rnd, breaks=20, plot=F)

d_dat <- data.frame(x=d_histinfo$mids, y=d_histinfo$density/max(d_histinfo$density))
w_dat <- data.frame(x=w_histinfo$mids, y=w_histinfo$density/max(w_histinfo$density))
m_dat <- data.frame(x=m_histinfo$mids, y=m_histinfo$density/max(m_histinfo$density))
l_dat <- data.frame(x=l_histinfo$mids, y=l_histinfo$density/max(l_histinfo$density))

ggplot() + xlab("log-returns") + ylab("relative frequency") + #xlim(-1,1) +
  geom_line(data=d_dat, aes(x, y), color='red') +
  geom_line(data=w_dat, aes(x, y), color='blue') +
  geom_line(data=m_dat, aes(x, y), color='yellow') +
  #geom_line(data=l_dat, aes(x, y), color='green') +
  scale_y_log10()

# Rescaling
d_dat <- data.frame(x=d_histinfo$mids, y=d_histinfo$density/max(d_histinfo$density))
w_dat <- data.frame(x=w_histinfo$mids/(7^0.546), y=w_histinfo$density/max(w_histinfo$density))
m_dat <- data.frame(x=m_histinfo$mids/(30^0.546), y=m_histinfo$density/max(m_histinfo$density))

ggplot() + xlab("log-returns") + ylab("relative frequency") + #xlim(-1,1) +
  geom_line(data=d_dat, aes(x, y), color='red') +
  geom_line(data=w_dat, aes(x, y), color='blue') +
  geom_line(data=m_dat, aes(x, y), color='yellow') +
  #geom_line(data=l_dat, aes(x, y), color='green') +
  scale_y_log10()

#
# Section # Fractals ###################################
#
# https://www.stat.washington.edu/sites/default/files/files/reports/2010/tr577.pdf
# https://cran.r-project.org/web/packages/fractaldim/fractaldim.pdf

fd_dcl <- fd.estimate(as.numeric(daily_cl$NEM),
                  methods = list(list(name="variation", p.index=0.5),
                                 "variogram", "hallwood", "boxcount"),
                  window.size = length(daily_cl$NEM), plot.loglog = T, nlags = 10)

fd_wcl <- fd.estimate(as.numeric(weekly_cl$NEM),
                     methods = list(list(name="variation", p.index=0.5),
                                    "variogram", "hallwood", "boxcount"),
                     window.size = length(weekly_cl$NEM), plot.loglog = T, nlags = 10)

fd_mcl <- fd.estimate(as.numeric(monthly_cl$NEM),
                      methods = list(list(name="variation", p.index=0.5),
                                     "variogram", "hallwood", "boxcount"),
                      window.size = length(monthly_cl$NEM), plot.loglog = T, nlags = 10)


fd_lr <- fd.estimate(as.numeric(daily_lr$NEM),
                  methods = list(list(name="variation", p.index=0.5),
                                 "variogram", "hallwood", "boxcount"),
                  window.size = length(daily_cl$NEM), plot.loglog = T, nlags = 10)

fd_wlr <- fd.estimate(as.numeric(weekly_lr$NEM),
                      methods = list(list(name="variation", p.index=0.5),
                                     "variogram", "hallwood", "boxcount"),
                      window.size = length(weekly_cl$NEM), plot.loglog = T, nlags = 10)

fd_mlr <- fd.estimate(as.numeric(monthly_lr$NEM),
                      methods = list(list(name="variation", p.index=0.5),
                                     "variogram", "hallwood", "boxcount"),
                      window.size = length(monthly_lr$NEM), plot.loglog = T, nlags = 10)


# Section # Links #########################################################
# https://otexts.org/fpp2/index.html