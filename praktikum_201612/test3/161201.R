

library("platexpress")
source("~/programs/platexpress/R/platereader.R")
source("~/programs/platexpress/R/interfaces.R")

experiment <- "161130_Praktikum_RAJ11_Test_3"

setwd("~/work/CoilProject/experiments/plategrowth/praktikum_201612/test3/")

plate.file <- "IPTG_Testreihe_3.csv"
data.file <- c("20161201_20161201 Praktikum - pRAJ11  1_Absorbance.CSV",
               "20161201_20161201 Praktikum - pRAJ11  1_Fluorescence.CSV")

## 1A) PARSE PLATE LAYOUT FILE
## SANITY CHECK PLATE LAYOUT: open file in spread-sheet program;
## check separator
plate <- readPlateMap(file=plate.file,sep=";",fsep=";",
                      blank.id="LB", fields=c("strain","IPTG","aTc"))

## 1B) PARSE DATA
## SANITY CHECK DATA: open file in spread-sheet program; check time
## format (see R strptime); and how many lines to skip!
raw <- readPlateData(file=data.file, type="BMG",time.conversion=1/3600)

groups <- getGroups(plate,by=c("strain"))
groups2 <- getGroups(plate,by=c("strain","IPTG","aTc"))
## TODO: allow better group sorting?
groups3 <- getGroups(plate,by=c("IPTG","aTc"))

raw <- prettyData(raw,dids=c(OD="584",YFP="485/Em520"),
                  col=c("#000000","#004AFF"))
data <- correctBlanks(raw,plate)

viewGroups(data, groups=groups)
viewGroups(data, groups=groups, groups2=groups2,show.ci=F,lwd.orig=1)
viewGroups(raw, groups=groups3,groups2=groups2,show.ci=F,lwd.orig=1)
