
library("platexpress")

setwd("~/work/CoilProject/experiments/plategrowth/ecoli_ts_20161014")

plate <-readPlateMap("20161014_platemap.csv",fsep=";",
                     fields=c("strain","IPTG","blank"))

files <- c("20161014_20161014 IPTG mVenus Injection  1_Absorbance.CSV",
           "20161014_20161014 IPTG mVenus Injection  1_Fluorescence.CSV")
raw <- readPlateData(files,type="BMG",time.conversion=1/60)

## time in minutes
raw <- prettyData(raw, dids=c(OD="584",mVenus="485/Em520"),
                  colors=c("#333300",wavelength2RGB(600)))

## find a nice color
showSpectrum()
findWavelength(3)

## view data
vp <- viewPlate(raw)
vp <- viewPlate(raw,xlim=c(0,1500),xscale=TRUE)

## GET A SINGLE DATASET
od <- getData(raw,"OD",type="orig")
TIME <- raw[["OD"]]$time/60
Xt <- od[,"A8"]

plot(TIME, Xt)

## cut to growth range
rng <- TIME<1500
Xt <- Xt[rng]
TIME <- TIME[rng]

## look at data
plot(TIME, Xt)
## log it?
plot(TIME, log(Xt))

## cut to linear range of growth
rng <- TIME>600 & TIME < 900
xt <- Xt[rng]
time <- TIME[rng]

## look again at data
plot(time,xt)
plot(time,log(xt))

## DO LINEAR REGRESSION
## ln(X(t)) = mu * t + ln(X(0))
lfit <- lm(log(xt) ~ time)

## check quality of fit
summary(lfit)

## get parameters
x0.1 <- exp(coefficients(lfit)[1])
mu.1 <- coefficients(lfit)[2]

## plot
lines(time, mu.1*time + log(x0.1), col="red")

## DO NON-LINEAR REGRESSION, using
## the results of the linear fit as initial parameter guesses
data <- data.frame(time=time, xt=xt)
start <- list(mu=mu.1,x0=x0.1)
nlfit <- nls(xt ~ x0*exp(mu*time),data=data,start=start)

## check quality of fit
summary(nlfit)

## get parameters
mu.2 <- coefficients(nlfit)[1]
x0.2 <- coefficients(nlfit)[2]

## plot results
plot(TIME, x0.2 * exp(TIME*mu.2),col="green", lty=2)
points(TIME,Xt)
lines(TIME, x0.1 * exp(TIME*mu.1),col="red")

## CONCLUSION??

## lag-phase, stationary phase
## Monod equation : mechanistic, dynamic model
## logistic, Gompertz, etc: phenomenological model
## allows to retrieve lag, mu and capacity

## do it all in batch by growth
## plot growth rate against conditions
