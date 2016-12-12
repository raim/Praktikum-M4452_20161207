library(platexpress)
library(grofit)
library(growthcurver)

dpath <- "./praktikum_201612/20161207/"
setwd(dpath)

use.bmg <- FALSE ## select BMG or Synergy 

## read plate layout
plate <- readPlateMap("IPTG_Testreihe_3.csv",sep = ";",fsep = ";",blank.id = "LB",fields = c("strain","IPTG","aTc"))

groups <- getGroups(plate, c("strain"), verb=T) # SET TO TRUE!
groups2 <- getGroups(plate, c("strain","IPTG"), verb=T)
groups3 <- getGroups(plate, c("strain","aTc","IPTG"), verb=T)
groups4 <- getGroups(plate, c("strain","aTc"), verb=T)


## read BMG data
files <- c("BMG/20161206_20161201 Praktikum - pRAJ11  1_Absorbance.CSV","BMG/20161206_20161201 Praktikum - pRAJ11  1_Fluorescence.CSV") #alter Plattenleser
raw <- readPlateData(files, type="BMG", time.conversion=1/3600)
vp<-viewPlate(raw)
raw<-prettyData(raw,dids=c(OD="584",GFP="485/Em520"),colors=c("#000000",wavelength2RGB(600)))
names(raw)
raw <- correctBlanks(raw, plate)
#raw <- adjustBase(raw,base=0, each=FALSE)

## read Synergy data
files2 <- c("Synergy/161207_Praktikum_realtime_GFP.csv") #neuer Plattenleser
raw2 <- readPlateData(files2, type="Synergy",skip=44,sep=";",time.format="%H:%M:%S",time.conversion=1/3600)
raw2<-prettyData(raw2,dids=c(OD="600",GFP="GFP_50:480,520"),colors=c("#000000",wavelength2RGB(600)))
names(raw2)
raw2 <- correctBlanks(raw2, plate)
#raw2 <- adjustBase(raw2,base=0,each=FALSE)


## select data set
if ( use.bmg ) {
    data <- raw
    od.rng <- 3
    gfpod.ylim <- c(0,2e3)
} else {
    data <- raw2
    od.rng <- .6
    gfpod.ylim <- c(0,1e4)
}


## normalize by OD
fl <- getData(data,"GFP")
od <- getData(data,"OD")
data <- addData(data, ID="GFP/OD", dat=fl/od, col="#0095FF")


## interpolate to OD
od.data <- interpolatePlateData(cutData(data,c(0,15)),"OD")
flod <- getData(od.data,"GFP/OD")

## ... and calculate fold-change to ininduced
uninduced <- rowMeans(flod[, groups3[["Z1_RAJ11_0a_0I"]]],na.rm=TRUE)#??????
od.data <- addData(od.data, ID="GFP/OD/uninduced", dat=flod/uninduced, col="#AAAA00")


## plot results
pdf(paste("Results",ifelse(use.bmg,"_BMG","_Synergy"),".pdf",sep=""),
    width=10,height=5)
viewGroups(data, groups=groups, groups2=groups3, show.ci=T,lwd.orig=0, verb=F,embed=FALSE,ylims=list("GFP/OD"=gfpod.ylim))
viewGroups(od.data, groups=groups, groups2=groups3,dids=c("GFP","GFP/OD"), show.ci=T,lwd.orig=0, verb=F,embed=FALSE)
viewGroups(od.data, groups=groups, groups2=groups3,dids="GFP/OD/uninduced", show.ci=T,lwd.orig=0, verb=F,embed=FALSE)
par(mai=c(2,.5,0,0))
results <- boxData(od.data, did="GFP/OD/uninduced", rng=od.rng, groups=groups3, type="bar")
dev.off()

