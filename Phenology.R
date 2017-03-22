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
## ADD GDD data
dd= degree.days.mat(cbind(clim$Min, clim$Max),LDT=12)
inds= which(!is.na(clim[,"Min"])& !is.na(clim[,"Max"]) )
clim$dd=NA
clim$dd[inds]= apply( clim[inds,c("Min","Max")], MARGIN=1, FUN=degree.days.mat, LDT=12 )  

#restrict to March 1 to Aug 31, ordinal 60 - 243
#code summer days
clim$summer=NA
clim[which(clim$Julian>59 & clim$Julian<244),"summer"]<-1
#set degree days outside summer to 0
clim["summer"==1,"dd"]=0

#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd = cumsum(dd)) 
#need to deal with NAs?

#---------------------
#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )

hop.b1= read.csv( "B1_1959-2015_eggDiapause.csv" )
hop.c1= read.csv( "C1_1959-2015_eggDiapause.csv" )
hop.cha= read.csv( "CHA_1959-2012_eggDiapause.csv" )
#add old a1 data
hop.a1= read.csv( "A1data_old.csv" )

#add site
hop.b1$site="B1"
hop.c1$site="C1"
hop.cha$site="CHA"

#combine
hop= rbind(hop.a1,hop.b1,hop.c1,hop.cha)

#add time period
hop$period="initial"
hop[which(hop$year>1960) ,"period"]="resurvey"

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
hop[which(hop$species=="Melanoplus dodgei"),"species"]="Melanoplus boulderensis"
hop[which(hop$species=="Melanoplus bivitattus"),"species"]="Melanoplus bivittatus"
hop[which(hop$species=="Cratypedes neglectus"),"species"]="Cratypledes neglectus"

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
      mean = mean(Mean, na.rm=TRUE), sd = sd(Mean, na.rm=TRUE),cdd = max(cdd, na.rm=TRUE),
      sem = sd(Mean, na.rm=TRUE)/sqrt(length(na.omit(Mean))))
#better NA handling?

#temp
ggplot(data=clim1, aes(x=Year, y = mean, color=Site ))+geom_point()+geom_line() +theme_bw()
#sd
ggplot(data=clim1, aes(x=Year, y = sd, color=Site ))+geom_point()+geom_line() +theme_bw()
#cum dd
ggplot(data=clim1, aes(x=Year, y = cdd, color=Site ))+geom_point()+geom_line() +theme_bw()

#---------
# Jadult and GDDadult ~year by sites

#subset to dates with adults
hop1= hop[which(hop$in6>0),]

#subset to focal species
specs= c("Aeropedellus clavatus","Chloealtis abdominalis","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes")

hop1= hop1[which(hop1$species %in% specs ),]

#metrics across years
hop2= ddply(hop1, c("site", "year","species","period"), summarise,
             min = min(ordinal, na.rm=TRUE), mean = mean(ordinal, na.rm=TRUE) )

#Change years for plotting ### FIX
hop2[which(hop2$year==1958),"year"]= 1958+40
hop2[which(hop2$year==1959),"year"]= 1959+40
hop2[which(hop2$year==1960),"year"]= 1960+40

#min ordinal date
ggplot(data=hop2, aes(x=year, y = min, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date
ggplot(data=hop2, aes(x=year, y = mean, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()

#-------------
## calculate median across individuals
hop1= hop1[order(hop1$ordinal),]
#cumulative sum of individuals within groups
hop1 = hop1 %>% group_by(species,site,year) %>% arrange(species,site,year,ordinal) %>% mutate(csind = cumsum(in6))
#number of median individual
hop3 = hop1 %>% group_by(species,site,year) %>% arrange(species,site,year,ordinal) %>% mutate(medind = max(csind)/2)

#date of median individual
hop3$inddif= abs(hop3$medind-hop3$csind) #difference from median individual

hop4= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif),]))

#plot
ggplot(data=hop4, aes(x=year, y = ordinal, color=site ))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()

#-----------------
# match phenology to temp and dd
## refine dd metrics? restrict to summer?

clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
hop2$siteyear= paste(hop2$site, hop2$year, sep="")
hop4$siteyear= paste(hop4$site, hop4$year, sep="")

hop2$Tmean=NA
hop2$cdd=NA
match1= match(hop2$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop2$Tmean[matched]<- clim1$mean[match1[matched]]  
hop2$cdd[matched]<- clim1$cdd[match1[matched]]  

hop4$Tmean=NA
hop4$cdd=NA
match1= match(hop4$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop4$Tmean[matched]<- clim1$mean[match1[matched]]  
hop4$cdd[matched]<- clim1$cdd[match1[matched]]

#TEMP PLOTS
#min ordinal date by temp
hop2= hop2[order(hop2$Tmean),]
ggplot(data=hop2, aes(x=Tmean, y = min, color=site,shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date by temp
ggplot(data=hop2, aes(x=Tmean, y = mean, color=site,shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#median individual ordinal by temp
ggplot(data=hop4, aes(x=Tmean, y = ordinal, color=site,shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()

#------------------
#GDD PLOTS
#min ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = min, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = mean, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#median individual ordinal by dd
ggplot(data=hop4, aes(x=cdd, y = ordinal, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()

#====================================
#ADDITIONAL GDD PLOTS

#subset to years
clim2= clim[which(clim$Year %in% c(1958:1960,2006:2015)),]
#subset sites
clim2= clim2[which(clim2$Site %in% c("A1","B1","C1")),]

#add time period
clim2$period="initial"
clim2[which(clim2$Year>1960) ,"period"]="resurvey"

#plot GDD accumulation over time
ggplot(data=clim2, aes(x=Julian, y = cdd, color=as.factor(Year),linetype=period))+geom_line()+facet_grid(.~Site) +theme_bw()

#------------------------------
#plot initial, recent GDDs

hop2$min.cdd=NA
hop2$mean.cdd=NA

#gdd at 1st appearance
hop.m= paste(hop2$siteyear,hop2$min,sep="")
clim.m= paste(clim$Site, clim$Year, clim$Julian,sep="")
#clim2[which(clim$Site==hop2$site[1] & clim$Year==hop2$year[1] & clim$Julian==hop2$min[1]),]
match1= match(hop.m,clim.m)
matched= which(!is.na(match1))
hop2$min.cdd[matched]<- clim$cdd[match1[matched]]  

#gdd at peak
hop.m= paste(hop2$siteyear,round(hop2$mean),sep="")
match1= match(hop.m,clim.m)
matched= which(!is.na(match1))
hop2$mean.cdd[matched]<- clim$cdd[match1[matched]] 

#gdd at 1st appearance
ggplot(data=hop2, aes(x=year, y = min.cdd, color=species))+geom_point()+geom_line()+facet_grid(.~site) +theme_bw()
#gdd at peak
ggplot(data=hop2, aes(x=year, y = mean.cdd, color=species))+geom_point()+geom_line()+facet_grid(.~site) +theme_bw()
