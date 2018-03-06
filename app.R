# Antonis Pishias
# Financial Data and Statistics
# Assignment 2

# Section 1 # Importing ###################################

setwd("~/Sites/rstudio/fds2")
setwd("./data//cryptocurrencyCSV/")

# install.packages("xts")
# install.packages("devtools")
# library(devtools)
# install_github("cran/PerformanceAnalytics")
library(xts)
library(PerformanceAnalytics)

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

# Section 2 # Helper Functions ###################################

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

# Section 3 # Time-Horizons ###################################
# Daily ####################
daily_lr <- mergeAll(cryptos_daily_lr)
periodicity(daily_lr)

# Scatter Plots
dfc <- as.data.frame(daily_lr)
pairs(dfc,diag.panel = panel.hist, upper.panel = panel.cor)

# Time-series Plots
dfcTS <- ts(dfc, frequency = 365, start = c(2015,06,26))
# Pairs
plot(dfcTS, main= "Cryptocurrency Daily Log-Returns")


# Weekly ####################
weekly_lr <- mergeAll(cryptos_weekly_lr)
periodicity(weekly_lr)

# Scatter Plots
wfc <- as.data.frame(weekly_lr)
pairs(wfc,diag.panel = panel.hist, upper.panel = panel.cor)

# Time-series Plots
wfcTS <- ts(wfc, frequency = 52, start = c(2015,06,26))
# Pairs
plot(wfcTS, main= "Cryptocurrency Weekly Log-Returns")

# Monthly ####################
monthly_lr <- mergeAll(cryptos_monthly_lr)
periodicity(monthly_lr)

# Scatter Plots
mfc <- as.data.frame(monthly_lr)
pairs(mfc,diag.panel = panel.hist, upper.panel = panel.cor)

# Time-series Plots
mfcTS <- ts(mfc, frequency = 12, start = c(2015,06,26))
# Pairs
plot(mfcTS, main= "Cryptocurrency Monthly Log-Returns")

# Section 4 # ARIMA ###################################
# This section assumes normality
# 
# install.packages("astsa")
library(astsa)

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


# Section 5 # Export ###################################
daily_cl <- as.data.frame(log(mergeAll(cryptos_daily_close)))
weekly_cl <- as.data.frame(log(mergeAll(cryptos_weekly_close)))
monthly_cl <- as.data.frame(log(mergeAll(cryptos_monthly_close)))

setwd("./output/")
write.csv(daily_cl, "cryptos_daily_log_close.csv", row.names = TRUE)
write.csv(weekly_cl, "cryptos_weekly_log_close.csv", row.names = TRUE)
write.csv(monthly_cl, "cryptos_monthly_log_close.csv", row.names = TRUE)

write.csv(dfc, "cryptos_daily_lr.csv", row.names = TRUE)
write.csv(wfc, "cryptos_weekly_lr.csv", row.names = TRUE)
write.csv(mfc, "cryptos_monthly_lr.csv", row.names = TRUE)
setwd("./../")


