

library("platexpress")
source("~/programs/platexpress/R/platereader.R")
source("~/programs/platexpress/R/interfaces.R")

experiment <- "161130_Praktikum_RAJ11_Test_2"

setwd("~/work/CoilProject/experiments/plategrowth/praktikum_201612/test2/")

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
raw <- readPlateData(file=data.file, type="Synergy",skip=44,sep=";",
                     time.format="%H:%M:%S",time.conversion=1/3600)

groups <- getGroups(plate,by=c("strain"))
groups2 <- getGroups(plate,by=c("strain","IPTG","aTc"))
## TODO: allow better group sorting?
groups3 <- getGroups(plate,by=c("IPTG","aTc"))

raw <- prettyData(raw,dids=c(OD="600",YFP="GFP_50:480,520"),
                  col=c("#000000","#004AFF"))
data <- correctBlanks(raw,plate)

viewGroups(data, groups=groups)
viewGroups(data, groups=groups, groups2=groups2,show.ci=F,lwd.orig=1)
viewGroups(raw, groups=groups3,groups2=groups2,show.ci=F,lwd.orig=1)
