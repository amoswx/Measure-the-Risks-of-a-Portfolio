library(readxl)
library(tidyverse)
library(RND)
library(fOptions)
library(FMStable)
library(xts)
library(rmgarch)
#part1
df=read_excel("F567.2022.project.SPXoptions.xlsx",skip=4)
df=df[-21,]
df$date=c(77,77,77,77,77,77,77,77,77,119,119,119,62,62,49,49,49,49,49,49)
df$rate=c(0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.2,0.2,0.19,0.19,0.19,0.19,0.19,0.19)
for (i in 1:nrow(df)){
  df$`Implied volatility`[i]=GBSVolatility(df$Price[i],"c",3269.96,df$`Strike Price`[i],df$date[i]/365,df$rate[i],df$rate[i]-0.01639)
}
df=df%>%
  mutate(Delta=Quantity*Multiplier*GBSGreeks("Delta","c",3269.96,`Strike Price`,date/365,rate,rate-0.01639,`Implied volatility`),
         Gamma=Quantity*Multiplier*GBSGreeks("Gamma","c",3269.96,`Strike Price`,date/365,rate,rate-0.01639,`Implied volatility`),
         Theta=Quantity*Multiplier*GBSGreeks("Theta","c",3269.96,`Strike Price`,date/365,rate,rate-0.01639,`Implied volatility`),
         Vega=Quantity*Multiplier*GBSGreeks("Vega","c",3269.96,`Strike Price`,date/365,rate,rate-0.01639,`Implied volatility`))
write.table(df,file="table.csv",sep = ",")
S0=seq(0.8*3269.96,1.2*3269.96,0.4*3269.96/79)
value=rep(0,80)
df1=df
for (i in 1:80){
  
  for (j in 1:nrow(df)){
    value2=GBSOption("c",S0[i],df$`Strike Price`[j],df$date[j]/365,df$rate[j],0.01639,df$`Implied volatility`[j])
    df1$price[j]=df$Quantity[j]*df$Multiplier[j]*value2@price
  }
  value[i]=sum(df1$price)
}
plot(x=S0,y=value,type="l")
write.table(S0,file="table1.csv",sep = ",")
write.table(value,file="table2.csv",sep = ",")
#part2
SPX_history <- read.csv("SPX_History.csv")
VIX_history <- read.csv("VIX_History.csv")
SPX <- subset(SPX_history, select = c(Date,Close))
VIX <- subset(VIX_history, select = c(Date,Close))
colnames(SPX)[colnames(SPX) == "Date"] <- "LDate"
colnames(VIX)[colnames(VIX) == "Date"] <- "LDate"
SPX$Date <- as.Date(SPX$LDate,format="%m/%d/%Y")
SPX <- subset(SPX, select = c(Date,Close))
VIX$Date <- as.Date(VIX$LDate,format="%m/%d/%Y")
VIX <- subset(VIX, select = c(Date,Close))
colnames(SPX)[colnames(SPX) == "Close"] <- "SPX"
colnames(VIX)[colnames(VIX) == "Close"] <- "VIX"
SPXVIX <- merge(SPX, VIX, by="Date")
SPXVIX$SPXret[2:nrow(SPXVIX)] <- log(SPXVIX$SPX[2:nrow(SPXVIX)]/SPXVIX$SPX[1:(nrow(SPXVIX)-1)])
SPXVIX$VIXret[2:nrow(SPXVIX)] <- log(SPXVIX$VIX[2:nrow(SPXVIX)]/SPXVIX$VIX[1:(nrow(SPXVIX)-1)])
valuedate = which(SPXVIX$Date == as.Date("10/30/2020",format="%m/%d/%Y"))
returns <- SPXVIX[(valuedate-999):valuedate,c(1,4,5)]
uspec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                     distribution.model = "norm")
uspec2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                     distribution.model = "norm") 
fit.marg1 <- ugarchfit(spec = uspec1, data = returns[,2])
fit.marg2 <- ugarchfit(spec = uspec2, data = returns[,3])
coef(fit.marg1)
coef(fit.marg2)
sigma_marg1 <- sigma(fit.marg1)
marginspec <- multispec(c(uspec1, uspec2))
mspec <- dccspec(marginspec, dccOrder = c(1,1), model = "DCC", distribution = "mvnorm")
mod <- dccfit(mspec,returns[,2:3])
mod
coef(mod)
SPXcoef <- fit.marg1@fit$coef    #GARCH model parameter estimates
SPXsigma <- fit.marg1@fit$sigma  #standard deviations for each date
SPXz <- fit.marg1@fit$z          #GARCH shocks z for each date
VIXcoef <- fit.marg2@fit$coef    #GARCH model parameter estimates
VIXsigma <- fit.marg2@fit$sigma  #standard deviations for each date
VIXz <- fit.marg2@fit$z          #GARCH shocks z for each date
DCCcoef <- mod@mfit$coef        #DCC parameter estimates (ncludes the GARCH estimates)
DCCR <- mod@mfit$R              #DCC correlations for each date
DCCQ <- mod@mfit$Q              #DCC q's for each date
SPXmean <- 0
VIXmean <- VIXcoef[1]+VIXcoef[2]*returns[1000,3]
SPXvariance <- SPXcoef[1]+SPXcoef[2]*returns[1000,2]^2+SPXcoef[3]*SPXsigma[1000]^2
VIXvariance <- VIXcoef[3]+VIXcoef[4]*returns[1000,3]^2+VIXcoef[5]*VIXsigma[1000]^2  
rhobar <- (1/1000)*SPXz %*% VIXz
q11 <- 1.0 + DCCcoef[9]*(SPXz[1000]^2 - 1) + DCCcoef[10]*(DCCQ[[1000]][1,1]-1)
q22 <- 1.0 + DCCcoef[9]*(VIXz[1000]^2 - 1) + DCCcoef[10]*(DCCQ[[1000]][2,2]-1)
q12 <- rhobar + DCCcoef[9]*(SPXz[1000]*VIXz[1000] - rhobar) + DCCcoef[10]*(DCCQ[[1000]][1,2]-rhobar)
rho = q12/sqrt(q11*q22)
delta <- -1126522 #The delta and vega I calculated and solution have a little different.
vega <- -2270796820 
portexpectedchange <- vega*(SPXVIX[valuedate,3]/100)*VIXmean
portvariance = delta^2*SPXVIX[valuedate,2]^2*SPXvariance + vega^2*(SPXVIX[valuedate,3]/100)^2*VIXvariance + delta*vega*SPXVIX[valuedate,2]*(SPXVIX[valuedate,3]/100)*rho*sqrt(SPXvariance*VIXvariance)
portsd <- sqrt(portvariance)
DNVaR = -1.645 *(portexpectedchange - portsd)
#Part3
dfa=returns%>%
  mutate(SPXsigma=SPXsigma,SPXz=SPXz,VIXsigma=VIXsigma,
         VIXz=VIXz,DCCR=lapply(lapply(DCCR,head,c(2,1)),tail,1),
         DCCQ11=lapply(DCCQ,head,c(1,1)),
         DCCQ22=lapply(DCCQ,tail,c(1,1)),
         DCCQ12=lapply(lapply(DCCQ,head,c(2,1)),tail,1))
z=dfa[,c(5,7)]
omega1=SPXcoef[1]
alpha1=SPXcoef[2]
beta1=SPXcoef[3]
sigma1=omega1/(1-alpha1-beta1)
omega2=VIXcoef[3]
alpha2=VIXcoef[4]
beta2=VIXcoef[5]
sigma2=omega2/(1-alpha2-beta2)
alphac=DCCcoef[9]
betac=DCCcoef[10]
q12=sum(dfa$SPXz*dfa$VIXz)/nrow(dfa)
dfa=as.numeric(dfa)
VIXmean <- VIXcoef[1]+VIXcoef[2]*returns[1000,3]
optiondata <- read.csv("optiondata.csv")
Strike     = optiondata[,2]
tau        = optiondata[,10]
rate       = optiondata[,11]
volatility = optiondata[,12]
Quantity   = optiondata[,5]
Price      = optiondata[,6]
value      = optiondata[,8]
S <- rep(1,20)*SPXVIX[valuedate,2]
q <- rep(1,20)*0.0163917
b = rate - q
optionPrice <- GBSOption("c", S=S, X=Strike, Time=tau, r=rate, b=b, sigma=volatility, title = NULL, 
                         description = NULL)@price
portValue <- 100 * Quantity %*% optionPrice
sigma1=sqrt(omega1+alpha1*dfa$SPXret[1000]^2+beta1*dfa$SPXsigma[1000]^2)
sigma2=sqrt(omega2+alpha2*dfa$VIXret[1000]^2+beta2*dfa$VIXsigma[1000]^2)
return1=sigma1*SPXz
return2=VIXmean+sigma2*VIXz
newtau = tau - 3/365
newPortValue <- rep(0,1000)
PL <- rep(0,1000)
for(i in 1:1000){
  newS <- SPXVIX[valuedate,2] * exp(return1[i]) 
  newVolatility1 <- volatility *exp(return2[i])
  newPrice <- GBSOption("c", S=newS, X=Strike, Time=newtau, r=rate, b=b, sigma=newVolatility1, title = NULL, 
                        description = NULL)@price
  newPortValue1[i] <- 100 * Quantity %*% newPrice
  PL[i] <- newPortValue1[i]-portValue
}
sum(is.na(PL1))
FHSVaR = -quantile(PL, probs = 0.05)
#part4
dfa=returns%>%
  mutate(SPXsigma=SPXsigma,SPXz=SPXz,VIXsigma=VIXsigma,
         VIXz=VIXz,DCCR=lapply(lapply(DCCR,head,c(2,1)),tail,1),
         DCCQ11=lapply(DCCQ,head,c(1,1)),
         DCCQ22=lapply(DCCQ,tail,c(1,1)),
         DCCQ12=lapply(lapply(DCCQ,head,c(2,1)),tail,1))
z=dfa[,c(5,7)]
omega1=SPXcoef[1]
alpha1=SPXcoef[2]
beta1=SPXcoef[3]
sigma1=omega1/(1-alpha1-beta1)
omega2=VIXcoef[3]
alpha2=VIXcoef[4]
beta2=VIXcoef[5]
sigma2=omega2/(1-alpha2-beta2)
alphac=DCCcoef[9]
betac=DCCcoef[10]
q12=sum(dfa$SPXz*dfa$VIXz)/nrow(dfa)
dfa=as.numeric(dfa)
VIXmean <- VIXcoef[1]+VIXcoef[2]*returns[1000,3]
dfb=as.data.frame(matrix(nrow=21,ncol=1))
names(dfb)=c("sigma1")
dfc=as.data.frame(matrix(nrow=10000,ncol=1))
names(dfc)=c("return1")
for (j in 1:10000){
  for (i in 1:21){
    if (i == 1){
      dfb$sigma1[1]=sqrt(omega1+alpha1*dfa$SPXret[1000]^2+beta1*dfa$SPXsigma[1000]^2)
      dfb$sigma2[1]=sqrt(omega2+alpha2*dfa$VIXret[1000]^2+beta2*dfa$VIXsigma[1000]^2)
      dfb$return1[1]=dfb$sigma1[1]*sample(dfa$SPXz,1,replace=TRUE)
      dfb$return2[1]=VIXmean+dfb$sigma2[1]*sample(dfa$VIXz,1,replace=TRUE)
    }
    else{
      dfb$sigma1[i]=sqrt(omega1+alpha1*dfb$return1[i-1]^2+beta1*dfb$sigma1[i-1]^2)
      dfb$sigma2[i]=sqrt(omega2+alpha2*dfb$return2[i-1]^2+beta2*dfb$sigma2[i-1]^2)
      VIXmean=VIXcoef[1]+VIXcoef[2]*returns[i-1,3]
      dfb$return1[i]=dfb$sigma1[i]*sample(dfa$SPXz,1,replace=TRUE)
      dfb$return2[i]=VIXmean+dfb$sigma2[i]*sample(dfa$VIXz,1,replace=TRUE)
    }
    dfc$return1[j]=sum(dfb$return1)
    dfc$return2[j]=sum(dfb$return2)
  }
}
simulatedSPX <- SPXVIX[valuedate,2] * exp(dfc$return1) 
simulatedVIX <- SPXVIX[valuedate,3] * exp(dfc$return2)
newtau1 = tau - 31/365
newPortValue1 <- rep(0,10000)
PL1 <- rep(0,10000)
for(i in 1:10000){
  newS <- SPXVIX[valuedate,2] * exp(dfc$return1[i]) 
  newVolatility1 <- volatility *exp(dfc$return2[i])
  newPrice <- GBSOption("c", S=newS, X=Strike, Time=newtau, r=rate, b=b, sigma=newVolatility1, title = NULL, 
                        description = NULL)@price
  newPortValue1[i] <- 100 * Quantity %*% newPrice
  PL1[i] <- newPortValue1[i]-portValue
}
sum(is.na(PL1))
FHSVaR1 = -quantile(PL1, probs = 0.05)
FHSVaR1
