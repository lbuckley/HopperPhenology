#load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(zoo)

#SUMMARY OF DATA PROCESSING
#1.	Fixed reconstruction of B1 2009 data, previously had issue with column selection
#2.	Coarsely reconstructed A1 data for 2009 and 2010 using data from A1, C1, NOAA, and Gross Reservoir
#3. Linearly interpolated data with a maximum gap of 6: B1 2012, C1 2013, C1 2015, A1 2010 

#source degree days function
#setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
setwd("/Users/laurenbuckley/HopperPhenology/")
source("degreedays.R")

#--------------------------------------
#fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

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
setwd( paste(fdir, "climate/Recent/", sep="") )  

#columns: "Site","Date","Year","Month","Julian","Max","Min","Mean"

#A1
clim.a1.2011= read.csv( "A1_2011.csv" )
clim.a1.2011$Site="A1"
clim.a1.2011$Date=NA
clim= rbind(clim, clim.a1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#B1
clim.b1.2011= read.csv( "B1_2011.csv" )
clim.b1.2011$Site="B1"
clim.b1.2011$Date=NA
clim= rbind(clim, clim.b1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

clim.b1.2015= read.csv( "B1_2015.csv" )
#calculate Julian
tmp <- as.POSIXlt(clim.b1.2015$Date.Time..GMT.07.00, format = "%m/%d/%Y %H:%M")
clim.b1.2015$Julian= tmp$yday+1
#calc min and max
clim.b1.2015= clim.b1.2015 %>% group_by(Year,Julian) %>% summarise(Date=Date.Time..GMT.07.00[1],Min= min(Temp_C),Max= max(Temp_C) )
clim.b1.2015= as.data.frame(clim.b1.2015)
clim.b1.2015$Site="B1"
clim.b1.2015$Mean=NA
clim.b1.2015$Month=NA
clim= rbind(clim, clim.b1.2015[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

clim.b1.2009.tmax= read.csv( "B1_2009_Tmax.csv" )
clim.b1.2009.tmin= read.csv( "B1_2009_Tmin.csv" )
clim.b1.2010.tmax= read.csv( "B1_2010_Tmax.csv" )
clim.b1.2010.tmin= read.csv( "B1_2010_Tmin.csv" )
clim.b1.2012.tmax= read.csv( "B1_2012_Tmax.csv" )
clim.b1.2012.tmin= read.csv( "B1_2012_Tmin.csv" )
clim.b1.2013.tmax= read.csv( "B1_2012_Tmax.csv" )

clim.b1.2009.tmax$Min= clim.b1.2009.tmin$Min
clim.b1.2010.tmax$Min= clim.b1.2010.tmin$Min
clim.b1.2012.tmax$Min= clim.b1.2012.tmin$Min
clim.b1.2013.tmax$Min=NA

clim.b1.20122013= rbind(clim.b1.2009.tmax, clim.b1.2010.tmax, clim.b1.2012.tmax, clim.b1.2013.tmax)
clim.b1.20122013$Site="B1"
clim.b1.20122013$Date=NA
clim.b1.20122013$Mean=NA
clim= rbind(clim, clim.b1.20122013[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#C1
clim.c1.20092012= read.csv( "C1_2009_2012.csv" )
clim.c1.2014= read.csv( "C1_2014.csv" )

clim.c1.2012= subset(clim.c1.20092012, clim.c1.20092012$Year>2011)
clim.c1.2012$Site="C1"
clim= rbind(clim, clim.c1.2012[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

clim.c1.2014$Site="C1"
clim.c1.2014$Month=NA
clim= rbind(clim, clim.c1.2014[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#------
#load and process LTER data
setwd( paste(fdir, "climate/LTERdownload_Recent/", sep="") )  

#A1
a1= read.csv("a-1hobo.hourly.jm.data.csv") #2013-2014
#calc min and max
a1= a1 %>% group_by(year,jday) %>% summarise(Date=date[1],Min= min(airtemp_avg),Max= max(airtemp_avg) )
a1= as.data.frame(a1)
a1$Site="A1"
a1$Year=a1$year
a1$Julian=a1$jday
a1$Mean=NA
a1$Month=NA

clim= rbind(clim, a1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#B1
b1= read.csv("b-1hobo.hourly.jm.data.csv") #2013-2014
#calc min and max
b1= b1 %>% group_by(year,jday) %>% summarise(Date=date[1],Min= min(airtemp_avg),Max= max(airtemp_avg) )
b1= as.data.frame(b1)
b1$Site="B1"
b1$Year=b1$year
b1$Julian=b1$jday
b1$Mean=NA
b1$Month=NA

clim= rbind(clim, b1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#C1
c1= read.csv("c-1tdayv.ml.data.csv")
#calculate Julian
tmp <- as.POSIXlt(c1$Date, format = "%m/%d/%Y")
c1$Julian= tmp$yday
c1$Site="C1"
c1$Month=NA
#subset to missing years
c1= subset(c1, c1$Year %in% c(2013,2015))

clim= rbind(clim, c1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#D1
d1= read.csv("d-1cr23x-cr1000.daily.ml.data_up.csv")
d1$Site="D1"
d1$Month=NA
#subset to missing years
d1= subset(d1, d1$Year>2008)

clim= rbind(clim, d1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean")])

#=====================================
# Remove duplicates

clim$Site_Year_Julian= paste(clim$Site, clim$Year, clim$Julian, sep="")
# average across duplicates
clim1= clim %>% group_by(Site,Date,Year,Month,Julian) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )
clim1= as.data.frame(clim1)

#replace NaN
clim1[clim1 == "NaN"] = "NA"
#change temps to numeric
clim1[,"Min"]= as.numeric(clim1[,"Min"])
clim1[,"Max"]= as.numeric(clim1[,"Max"])
clim1[,"Mean"]= as.numeric(clim1[,"Mean"])

#------------------------
#Plot

#summer means
clim2= clim1[which(clim1$Julian>59 & clim1$Julian<244),]
clim2 = clim2 %>% group_by(Year,Site) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )

ggplot(data=clim2, aes(x=Year, y = Min, color=Site ))+geom_line() +theme_bw()

#=====================================
# WRITE OUT DATA
setwd( paste(fdir, "climate", sep="") )   
#write.csv(clim1, "AlexanderClimateAll.csv")

#check counts of data
years= c(1958:1960, 2006:2015)

climsub=clim1[clim1$Year %in% years,]
#counts across sites years
climsub= climsub %>% group_by(Site,Year) %>% summarise(Min= length(na.omit(Min)),Max= length(na.omit(Max)),Mean= length(na.omit(Mean)) )
climsub= as.data.frame(climsub)

#-------------------------------
#CHECK AND INTERPOLATE CLIMATE DATA

setwd( paste(fdir, "climate", sep="") )   
clim1= read.csv("AlexanderClimateAll.csv")

# Subset climate data to sites and years for which we have grasshopper data
sites <- c("A1", "B1", "C1", "NOAA")
years <- c(1958:1960, 2006:2015)
allClim <- droplevels(clim1[clim1$Site %in% sites & clim1$Year %in% years,])

# Subset to ordinal dates relevant for grasshopper season (March 1 to Aug 31 = ordinal 60 to 243) #extend to 290
allClim <- allClim[allClim$Julian > 59 & allClim$Julian < 290,]

#set up data frame with all combinations
clim1 = expand.grid(Site=sites, Julian = 60:290, Year = years)
#set up columns
clim1$Max=NA; clim1$Min=NA; clim1$Mean=NA
#columns for matching
clim1$sjy= paste(clim1$Site, clim1$Julian, clim1$Year, sep="_")
allClim$sjy= paste(allClim$Site, allClim$Julian, allClim$Year, sep="_")
#match
match1= match(clim1$sjy,allClim$sjy)
matched= which(!is.na(match1))

clim1[matched,c("Site","Julian","Year","Max","Min","Mean")]= allClim[match1[matched], c("Site","Julian","Year","Max","Min","Mean")]
# Re-order Climate data
clim1 <- clim1[order(clim1$Site, clim1$Year, clim1$Julian),]

#---------------------
#Update climate data
setwd( paste(fdir, "climate/Recent/", sep="") ) 

#Fix B1 2009 data
clim.b1.2009.tmax= read.csv( "B1_2009_Tmax.csv" )
clim.b1.2009.tmin= read.csv( "B1_2009_Tmin.csv" )

clim.b1.2009.tmax$Min= clim.b1.2009.tmin$Min

clim.b1.2009.tmax$Site="B1"
clim.b1.2009.tmax$sjy= paste(clim.b1.2009.tmax$Site, clim.b1.2009.tmax$Julian, clim.b1.2009.tmax$Year, sep="_")
#match
match1= match(clim1$sjy, clim.b1.2009.tmax$sjy)
matched= which(!is.na(match1))

clim1[matched,c("Max","Min")]= clim.b1.2009.tmax[match1[matched], c("Max","Min")]

#----
#ADD C1 2015
clim.c1.2015= read.csv( "C1_2015_CHECK.csv" )
#calculate Julian
tmp <- as.POSIXlt(clim.c1.2015$C1.Temp.data, format = "%m/%d/%Y %H:%M")
clim.c1.2015$Julian= tmp$yday+1
#calc min and max
clim.c1.2015= clim.c1.2015 %>% group_by(Julian) %>% summarise(Min= min(Tc),Max= max(Tc) )
clim.c1.2015= as.data.frame(clim.c1.2015)
clim.c1.2015$Site="C1"
clim.c1.2015$Year=2015
clim.c1.2015$sjy= paste(clim.c1.2015$Site, clim.c1.2015$Julian, clim.c1.2015$Year, sep="_")
#match
match1= match(clim1$sjy, clim.c1.2015$sjy)
matched= which(!is.na(match1))

clim1[matched,c("Max","Min")]= clim.c1.2015[match1[matched], c("Max","Min")]

#----
#ADD RECONSTRUCTED A1 data
clim.a1= read.csv("A1_2009_2010_reconstruct.csv")

clim.a1$Site="A1"
clim.a1$sjy= paste(clim.a1$Site, clim.a1$Julian, clim.a1$Year, sep="_")
#match
match1= match(clim1$sjy, clim.a1$sjy)
matched= which(!is.na(match1))

clim1[matched,c("Max","Min")]= clim.a1[match1[matched], c("Max","Min")]

#------------------------

clim.nas= clim1[(is.na(clim1$Min) | is.na(clim1$Max)),] 
clim.nas$Year= as.factor(clim.nas$Year)
#counts by sites, years
clim.nas2= clim.nas %>% group_by(Site, Year) %>% summarise(count=length(Julian) )

# Interpolate data
inds= which(clim1$Site=="B1"&clim1$Year=="2012")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=5)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=5)
#approx first value
clim1$Min[inds[1]]= clim1$Min[inds[2]]
clim1$Max[inds[1]]= clim1$Max[inds[2]]
#-------
inds= which(clim1$Site=="C1"&clim1$Year=="2013")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=5)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=5)
#--------
inds= which(clim1$Site=="C1"&clim1$Year=="2015")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=6)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=6)
#--------
inds= which(clim1$Site=="A1"&clim1$Year=="2010")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=6)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=6)

#--------------------
#Add degree days
clim= clim1

#Use NOAA for CHA
clim$Site= as.character(clim$Site)
clim[which(clim$Site=="NOAA"),"Site"]<-"CHA"
clim$Site= as.factor(clim$Site)

#--------------------------------
## ADD GDD data
dd= degree.days.mat(cbind(clim$Min, clim$Max),LDT=12)
inds= which(!is.na(clim[,"Min"])& !is.na(clim[,"Max"]) )
clim$dd=NA
clim$dd[inds]= apply( clim[inds,c("Min","Max")], MARGIN=1, FUN=degree.days.mat, LDT=12 )  

#-----------
# Calculate additional DD metrics
#Summer 60-243
#June 152-181
#July 182-212
#Aug 213-243
#Early 60-151
#Mid 60-181

clim$dd_sum=clim$dd
clim[which(clim$Julian<60 | clim$Julian>243),"dd_sum"]<-0

#add fall
clim$dd_sumfall=clim$dd
clim[which(clim$Julian<60),"dd_sumfall"]<-0

clim$dd_june=clim$dd
clim[which(clim$Julian<152 | clim$Julian>181),"dd_june"]<-0

clim$dd_july=clim$dd
clim[which(clim$Julian<182 | clim$Julian>212),"dd_july"]<-0

clim$dd_aug=clim$dd
clim[which(clim$Julian<213 | clim$Julian>243),"dd_aug"]<-0

clim$dd_early=clim$dd
clim[which(clim$Julian<60 | clim$Julian>151),"dd_early"]<-0

clim$dd_mid=clim$dd
clim[which(clim$Julian<60 | clim$Julian>181),"dd_mid"]<-0

#---------------------
#species specific cdd: month before mean ordinal date
clim$dd_ac=clim$dd
clim[which(clim$Julian<(195-30) | clim$Julian>(195+30) ),"dd_ac"]<-0
clim$dd_mb=clim$dd
clim[which(clim$Julian<(208-30) | clim$Julian>(208+30) ),"dd_mb"]<-0
clim$dd_ms=clim$dd
clim[which(clim$Julian<(224-30) | clim$Julian>(224+30) ),"dd_ms"]<-0

#=========================
#WRITE OUT
setwd( paste(fdir, "climate", sep="") )   
write.csv(clim,"AlexanderClimateAll_filled_May2022.csv")

#--------------------------
#CHECK DATA

#B1
clim.b1= clim[clim$Site=="B1" & clim$Julian %in% 60:243,]
clim.b1$Year= as.factor(clim.b1$Year)

ggplot(data=clim.b1, aes(x=Julian, y = Min, color=Year))+geom_smooth() +theme_bw()
ggplot(data=clim.b1, aes(x=Julian, y = Max, color=Year))+geom_smooth() +theme_bw()

#------
#A1
clim=clim1
clim.a1= clim[clim$Site=="A1" & clim$Julian %in% 60:243,]
clim.a1$Year= as.factor(clim.a1$Year)

ggplot(data=clim.a1, aes(x=Julian, y = Min, color=Year))+geom_smooth() +theme_bw()
ggplot(data=clim.a1, aes(x=Julian, y = Max, color=Year))+geom_smooth() +theme_bw()

