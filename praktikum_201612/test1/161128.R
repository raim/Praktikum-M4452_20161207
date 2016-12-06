library(platexpress)

experiment <- "161128_Praktikum_RAJ11_Test_1"

setwd("~/work/hhu_2015/uebung_201612/Praktikum-M4452_20161207/praktikum_201612/test1/")

plate.file <- paste(experiment,"_layout.csv",sep="")
data.file <- paste(experiment,".csv",sep="")

## 1A) PARSE PLATE LAYOUT FILE
## SANITY CHECK PLATE LAYOUT: open file in spread-sheet program;
## check separator
plate <- readPlateMap(file=plate.file,sep="\t",fsep=";",
                      blank.id="LB", fields=c("strain","IPTG","aTc"))

## 1B) PARSE DATA
## SANITY CHECK DATA: open file in spread-sheet program; check time
## format (see R strptime); and how many lines to skip!
raw <- readPlateData(file=data.file, type="Synergy",skip=55,sep=";",
                     time.format="%H:%M:%S",time.conversion=1/3600)


## CHECK DATA
viewPlate(raw)

## -> all wells seem to be Ok, best wavelength is YFP_50:500,535

raw <- prettyData(raw, dids=c(OD="600", YFP="YFP_50:500,535"))

## correct blanks
data <- correctBlanks(raw, plate)

## get Groupings
g1 <- getGroups(plate,"strain")
g2 <- getGroups(plate,c("strain","IPTG"))
g3 <- getGroups(plate,c("strain","IPTG","aTc"))

viewGroups(data,groups=g1)
viewGroups(data,groups=g1,groups2=g2)
viewGroups(data,groups=g2,groups2=g3)
viewGroups(data,groups=g1,groups2=g3)
