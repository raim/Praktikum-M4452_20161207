library(platexpress)
vignette("platexpress")
setwd('/home/tim/quabi/synthetische_biologie')


plate <- readPlateMap("Praktikum-M4452_20161207/praktikum_201612/test3/IPTG_Testreihe_3.csv",sep = ";", fsep = ";", blank.id="LB",fields=c("strain","IPTG","ATC"))

# Gruppierungen der Platte
groups <- getGroups(plate, c("strain"),verb=F)
groups2 <- getGroups(plate,c("strain", "IPTG"),verb=F)
groups3 <- getGroups(plate,c("strain", "IPTG","ATC"),verb=F)
groups4 <- getGroups(plate,c("ATC"),verb=F)


#Daten aus Platereader Rainer
files1 <- c("Praktikum-M4452_20161207/praktikum_201612/20161207_BMG/20161206_20161201 Praktikum - pRAJ11  1_Absorbance.CSV","Praktikum-M4452_20161207/praktikum_201612/20161207_BMG/20161206_20161201 Praktikum - pRAJ11  1_Fluorescence.CSV") 

raw1 <- readPlateData(files1,type="BMG", time.conversion=1/60)
# Rohansicht der Daten
vp1 = viewPlate(raw1)

raw1 <- prettyData(raw1,dids=c(OD="584",mVenus="485/Em520"))
viewGroups(raw1,groups=groups, groups2=groups2,verb=F)
viewGroups(raw1,groups=groups, groups2=groups3,verb=F)
#Daten aus Platereader Kollmann
files2 <- c("161207_Praktikum_realtime_GFP.csv")
raw2 <- readPlateData(files2,type="Synergy", skip=44, time.conversion=1/60, time.format="%H:%M:%S")

vp2 <- viewPlate(raw2)
names(raw2)
raw2 <- prettyData(raw2,dids=c(OD="600",mVenus="GFP_50:480,520"))

viewGroups(raw2,groups=groups, groups2=groups2,verb=F)
viewGroups(raw2,groups=groups, groups2=groups3,verb=F)
viewGroups(raw2, groups2=groups4,verb=F)

# Prepare Data (blanks, cuts, etc. )
# Weiterarbeiten mit den Sympy Daten

raw3 <- correctBlanks(raw2,plate)

# Berechnen der Wachstumsraten (ein Beispiel)
od = getData(raw3,"OD")
TIME <- raw3$Time
Xt <- od[,"C7"]
plot(TIME,Xt)
plot(TIME,log(Xt))

rng <- TIME<200
xt <- Xt[rng]
t = TIME[rng]
par(mfcol=c(1,2))
plot(t,xt)
plot(t,log(xt))

# lm model und nl model

lmfit = lm(log(xt) ~ t)

x0.1 = exp(coefficients(lmfit)[1])
mu.1 = coefficients(lmfit)[2]

plot(t,xt)
lines(t,x0.1*exp(mu.1*t))

dat =data.frame(time=t,xt=xt)
start=list(mu=mu.1,x0=x0.1)
nlfit <- nls(xt ~ x0*exp(mu*time),data=dat,start=start)

plot(t,xt)
lines(t,coefficients(nlfit)['x0']*exp(coefficients(nlfit)['mu']*t))




