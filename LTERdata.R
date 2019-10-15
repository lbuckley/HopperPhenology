library(dplyr)
library(zoo)

fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"
setwd("/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/ForLTER/final/")

# #source degree days function
# setwd("/Users/laurenbuckley/HopperPhenology/")
# source("degreedays.R")

#--------------------------------------
#HOPPERS

hop.dat= read.csv("hoppers.cn.data_predate.csv")
hop.dat$date= as.character(format( as.Date(hop.dat$date, format="%m/%d/%y"), format="%Y-%m-%d") )
#Fix 2059... dates
hop.dat$date= gsub("2058","1958",hop.dat$date)
hop.dat$date= gsub("2059","1959",hop.dat$date)
hop.dat$date= gsub("2060","1960",hop.dat$date)

write.csv(hop.dat, "hoppers.cn.data.csv")

#--------------------------------------
#CLIMATE

#SUMMARY OF DATA PROCESSING
#1.	Fixed reconstruction of B1 2009 data, previously had issue with column selection
#2.	Coarsely reconstructed A1 data for 2009 and 2010 using data from A1, C1, NOAA, and Gross Reservoir
#3. Linearly interpolated data with a maximum gap of 6: B1 2012, C1 2013, C1 2015, A1 2010 

#load climate data
setwd("/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/climate/")

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

#combine
clim= rbind(clim.a1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","FlagMax","FlagMin")], clim.b1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","FlagMax","FlagMin")], clim.c1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","FlagMax","FlagMin")], clim.d1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","FlagMax","FlagMin")], clim.noaa[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","FlagMax","FlagMin")])
clim$Source= "measured"
#account for flags
inds= which(!clim$FlagMin=="" | !clim$FlagMax=="")
clim$Source[inds]="modelled"
#drop flag columns
clim= clim[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")]

#------
#load recent climate data from Cesar
setwd("/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/climate/Recent/")

#columns: "Site","Date","Year","Month","Julian","Max","Min","Mean"

#A1
clim.a1.2011= read.csv( "A1_2011.csv" )
clim.a1.2011$Site="A1"
clim.a1.2011$Date=NA
clim.a1.2011$Source= "measured"
clim= rbind(clim, clim.a1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

#B1
clim.b1.2011= read.csv( "B1_2011.csv" )
clim.b1.2011$Site="B1"
clim.b1.2011$Date=NA
clim.b1.2011$Source= "measured"
clim= rbind(clim, clim.b1.2011[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

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
clim.b1.2015$Source= "measured"
clim= rbind(clim, clim.b1.2015[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

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
clim.b1.20122013$Source= "measured"
clim= rbind(clim, clim.b1.20122013[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

#C1
clim.c1.20092012= read.csv( "C1_2009_2012.csv" )
clim.c1.2014= read.csv( "C1_2014.csv" )

clim.c1.2012= subset(clim.c1.20092012, clim.c1.20092012$Year>2011)
clim.c1.2012$Site="C1"
clim.c1.2012$Source= "measured"
clim= rbind(clim, clim.c1.2012[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

clim.c1.2014$Site="C1"
clim.c1.2014$Month=NA
clim.c1.2014$Source= "measured"
clim= rbind(clim, clim.c1.2014[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

#------
#load and process LTER data
setwd("/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/climate/LTERdownload_Recent/")

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
a1$Source= "measured"

clim= rbind(clim, a1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

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
b1$Source= "measured"

clim= rbind(clim, b1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

#C1
c1= read.csv("c-1tdayv.ml.data.csv")
#calculate Julian
tmp <- as.POSIXlt(c1$Date, format = "%m/%d/%Y")
c1$Julian= tmp$yday
c1$Site="C1"
c1$Month=NA
#subset to missing years
c1= subset(c1, c1$Year %in% c(2013,2015))
c1$Source= "measured"

clim= rbind(clim, c1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

#D1
d1= read.csv("d-1cr23x-cr1000.daily.ml.data_up.csv")
d1$Site="D1"
d1$Month=NA
#subset to missing years
d1= subset(d1, d1$Year>2008)
d1$Source= "measured"

clim= rbind(clim, d1[,c("Site","Date","Year","Month","Julian","Max","Min","Mean","Source")])

#=====================================
# Remove duplicates

clim$Site_Year_Julian= paste(clim$Site, clim$Year, clim$Julian, sep="")
# average across duplicates
clim1= clim %>% group_by(Site,Date,Year,Month,Julian, Source) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )
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

# WRITE OUT DATA
#setwd( paste(fdir, "climate", sep="") )   
#write.csv(clim1, "AlexanderClimateAll.csv")

#=====================================
#CHECK AND INTERPOLATE CLIMATE DATA

#setwd( paste(fdir, "climate", sep="") )   
#clim1= read.csv("AlexanderClimateAll.csv")

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

clim1[matched,c("Site","Julian","Year","Max","Min","Mean","Source")]= allClim[match1[matched], c("Site","Julian","Year","Max","Min","Mean","Source")]
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
clim1[matched,"Source"]= "modelled"

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
clim1[matched,"Source"]= "measured"

#----
#ADD RECONSTRUCTED A1 data
clim.a1= read.csv("A1_2009_2010_reconstruct.csv")

clim.a1$Site="A1"
clim.a1$sjy= paste(clim.a1$Site, clim.a1$Julian, clim.a1$Year, sep="_")
#match
match1= match(clim1$sjy, clim.a1$sjy)
matched= which(!is.na(match1))

clim1[matched,c("Max","Min")]= clim.a1[match1[matched], c("Max","Min")]
clim1[matched,"Source"]= "modelled"

#------------------------

clim.nas= clim1[(is.na(clim1$Min) | is.na(clim1$Max)),] 
clim.nas$Year= as.factor(clim.nas$Year)
#counts by sites, years
clim.nas2= clim.nas %>% group_by(Site, Year) %>% summarise(count=length(Julian) )
clim.nas2= as.data.frame(clim.nas2)

# Interpolate data
inds= which(clim1$Site=="B1"&clim1$Year=="2012")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=5)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=5)
inds= which(clim1$Site=="B1"&clim1$Year=="2012"&is.na(clim1$Min) )
clim1$Source[inds] <- "modelled"
#approx first value
clim1$Min[inds[1]]= clim1$Min[inds[2]]
clim1$Max[inds[1]]= clim1$Max[inds[2]]
#-------
inds= which(clim1$Site=="C1"&clim1$Year=="2013")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=5)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=5)
inds= which(clim1$Site=="C1"&clim1$Year=="2013"&is.na(clim1$Min) )
clim1$Source[inds] <- "modelled"
#--------
inds= which(clim1$Site=="C1"&clim1$Year=="2015")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=6)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=6)
inds= which(clim1$Site=="C1"&clim1$Year=="2015"&is.na(clim1$Min) )
clim1$Source[inds] <- "modelled"
#--------
inds= which(clim1$Site=="A1"&clim1$Year=="2010")
clim1$Min[inds] <- na.approx(clim1$Min[inds], na.rm = F, maxgap=6)
clim1$Max[inds] <- na.approx(clim1$Max[inds], na.rm = F, maxgap=6)
inds= which(clim1$Site=="A1"&clim1$Year=="2010"&is.na(clim1$Min) )
clim1$Source[inds] <- "modelled"

#=========================
#WRITE OUT
setwd( paste(fdir, "climate", sep="") )   
#change Julian to Ordinal
names(clim1)[2]="Ordinal"
write.csv(clim1[,c(1:5,8)], "AlexanderClimateAll_filled_Oct2019.csv", row.names = FALSE)

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

