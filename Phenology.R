#load libraries
library(ggplot2)
library(plyr)

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
## ADD GDD data

#---------------------
#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )

hop.b1= read.csv( "B1_1959-2015_eggDiapause.csv" )
hop.c1= read.csv( "C1_1959-2015_eggDiapause.csv" )
hop.cha= read.csv( "CHA_1959-2012_eggDiapause.csv" )

#add site
hop.b1$site="B1"
hop.c1$site="C1"
hop.cha$site="CHA"

#combine
hop= rbind(hop.b1,hop.c1,hop.cha)

#--------------------------------------
#CHECK DATA

#check years
sort(unique(clim$Year)) #through 2014
sort(unique(hop$year)) #through 2015

#remove white space in species names
hop$species=trimws(hop$species, which = "both")
#check species
sort(unique(hop$species))
#Fix typos
hop[which(hop$species=="Melanoplus bouderensis"),"species"]="Melanoplus boulderensis"

#--------------------------------------
#ANALYSIS

#temperature plots
# temp, var, and GDD ~year by sites

#subset years to survey
clim1= clim[which(clim$Year %in% c(1959, 1960, 2006:2015))  ,]

#Change years for plotting ### FIX
clim1[which(clim1$Year==1959),"Year"]= 1959+40
clim1[which(clim1$Year==1960),"Year"]= 1960+40

#metrics across years
clim1= ddply(clim1, c("Site", "Year"), summarise,
      mean = mean(Mean, na.rm=TRUE), sd = sd(Mean, na.rm=TRUE),
      sem = sd(Mean, na.rm=TRUE)/sqrt(length(na.omit(Mean))))
#better NA handling?

#temp
ggplot(data=clim1, aes(x=Year, y = mean, color=Site ))+geom_point()+geom_line() +theme_bw()
#sd
ggplot(data=clim1, aes(x=Year, y = sd, color=Site ))+geom_point()+geom_line() +theme_bw()

#---------
# Jadult and GDDadult ~year by sites

#subset to dates with adults
hop1= hop[which(hop$in6>0),]

#subset to focal species
specs= c("Aeropedellus clavatus","Chloealtis abdominalis","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes")

hop1= hop1[which(hop1$species %in% specs ),]

#metrics across years
hop2= ddply(hop1, c("site", "year","species"), summarise,
             min = min(ordinal, na.rm=TRUE), mean = mean(ordinal, na.rm=TRUE) )

#Change years for plotting ### FIX
hop2[which(hop2$year==1959),"year"]= 1959+40
hop2[which(hop2$year==1960),"year"]= 1960+40

#min ordinal date
ggplot(data=hop2, aes(x=year, y = min, color=site ))+geom_point()+geom_line()+facet_wrap(~species, ncol=4) +theme_bw()
#mean ordinal date
ggplot(data=hop2, aes(x=year, y = mean, color=site ))+geom_point()+geom_line()+facet_wrap(~species, ncol=4) +theme_bw()

## calculate median across individuals

#-----------------
# match phenology to temp


