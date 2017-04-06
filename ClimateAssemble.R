#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

#--------------------------------------
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   

clim.a1= read.csv( "Climate_A1_1953-2014.csv" )
clim.b1= read.csv( "Climate_B1_1953-2014.csv" )
clim.c1= read.csv( "Climate_C1_1953-2014.csv" )
clim.d1= read.csv( "Climate_D1_1953-2014.csv" )
clim.noaa= read.csv( "Climate_noaa_1953-2014.csv" )

#add site
clim.a1$Site="A1"
clim.b1$Site="B1"
clim.c1$Site="C1"
clim.d1$Site="D1"
clim.noaa$Site="NOAA"

#combine ##currently omits flags
clim= rbind(clim.a1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.b1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.c1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.d1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")], clim.noaa[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#------
#load recent climate data from Cesar
setwd( paste(fdir, "climate\\Recent\\", sep="") )  

#columns: "Site","Date","Year","Month","Julian","Max","Min","Mean"

#A1
clim.a1.2011= read.csv( "A1_2011.csv" )
clim.a1.2011$Site="A1"
clim.a1.2011$Date=NA
clim= rbind(clim, clim.a1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#B1
clim.b1.2011= read.csv( "B1_2011.csv" )

clim.b1.2012.tmax= read.csv( "B1_2012_Tmax.csv" )
clim.b1.2012.tmin= read.csv( "B1_2012_Tmin.csv" )
clim.b1.2013.tmax= read.csv( "B1_2012_Tmax.csv" )

clim.b1.2012.tmax$Min= clim.b1.2012.tmin$Min
clim.b1.2013.tmax$Min=NA
clim.b1.20122013= rbind(clim.b1.2012.tmax,clim.b1.2013.tmax)

clim.b1.20122013
