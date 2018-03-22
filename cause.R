setwd("~/Sites/rstudio/fds2")

cryto_lr <- read.csv(".//data//output//all_crypto_log_returns.csv")
cryto_lr <- cryto_lr[-1]

library(aTSA)
stationary.test(cryto_lr$Bitcoin, method="kpss")
stationary.test(cryto_lr$NEM, method="kpss")
stationary.test(cryto_lr$Stellar, method="kpss")
stationary.test(cryto_lr$Litecoin, method="kpss")
stationary.test(cryto_lr$Siacon, method="kpss")

### Cause
library(lmtest)
# NEM causes Bitcoin????
grangertest(Bitcoin ~ NEM, order=3, data=cryto_lr)
# Bitcoin causes NEM????
grangertest(NEM ~ Bitcoin, order=3, data=cryto_lr)

grangertest(Stellar ~ NEM, order=3, data=cryto_lr)
grangertest(NEM ~ Stellar, order=3, data=cryto_lr)

grangertest(Litecoin ~ NEM, order=3, data=cryto_lr)
grangertest(NEM ~ Litecoin, order=3, data=cryto_lr)

grangertest(Siacon ~ NEM, order=3, data=cryto_lr)
grangertest(NEM ~ Siacon, order=3, data=cryto_lr)





grangertest(Ripple ~ NEM, order=3, data=cryto_lr)
grangertest(NEM ~ Ripple, order=3, data=cryto_lr)

grangertest(Monero ~ NEM, order=3, data=cryto_lr)
grangertest(NEM ~ Monero, order=3, data=cryto_lr)

grangertest(Dash ~ NEM, order=3, data=cryto_lr)
grangertest(NEM ~ Dash, order=3, data=cryto_lr)



grangertest(Ripple ~ Stellar, order=3, data=cryto_lr)
grangertest(Stellar ~ Ripple, order=3, data=cryto_lr)

grangertest(Bitcoin ~ Litecoin, order=3, data=cryto_lr)
grangertest(Litecoin ~ Bitcoin, order=3, data=cryto_lr)

# gc <- grangertest(Bitcoin ~ NEM, order=i, data=cryto_lr)
# 
print("Does Litecoin Grangerly Cause Bitcoin?")
for (i in 1:50) {
  gc <- grangertest(Bitcoin ~ Litecoin, order=i, data=cryto_lr)
  print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
}


print("Does Stellar Grangerly Cause Ripple?")
for (i in 1:50) {
  gc <- grangertest(Ripple ~ Stellar, order=i, data=cryto_lr)
  print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
}

print("Does Ripple Grangerly Cause Stellar?")
for (i in 1:50) {
  gc <- grangertest(Stellar ~ Ripple, order=i, data=cryto_lr)
  print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
}
# 
# print("Does Bitcoin Grangerly Cause NEM?")
# for (i in 1:50) {
#   gc <- grangertest(Bitcoin ~ NEM, order=i, data=cryto_lr)
#   print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
# }
# 
# print("Does Stellar Grangerly Cause NEM?")
# for (i in 1:50) {
#   gc <- grangertest(Stellar ~ NEM, order=i, data=cryto_lr)
#   print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
# }
# 
# print("Does Litecoin Grangerly Cause NEM?")
# for (i in 1:50) {
#   gc <- grangertest(Litecoin ~ NEM, order=i, data=cryto_lr)
#   print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
# }
# 
# print("Does Siacon Grangerly Cause NEM?")
# for (i in 1:50) {
#   gc <- grangertest(Siacon ~ NEM, order=i, data=cryto_lr)
#   print(paste(gc$Df[2],gc$`Pr(>F)`[2], gc$F[2], sep = "--"))
# }
