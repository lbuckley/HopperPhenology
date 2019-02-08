#Abbreviate to focus on summer GDDs and quantiles

#load libraries
library(ggplot2)
library(plyr) 
library(dplyr)
library(colorRamps)
require(cowplot)

sites= c("CHA", "A1", "B1", "C1", "D1")  #Redfox: 1574
elevs= c(1752, 2195, 2591, 3048, 3739)

#source degree days function
setwd("/Users/laurenbuckley/HopperPhenology/")
source("degreedays.R")

#======================================================
#READ DATA

fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled_May2018.csv")

#---------------------
#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd_sum = cumsum(dd_sum), cdd_sumfall = cumsum(dd_sumfall), cdd_july = cumsum(dd_july)) 
# cdd_june = cumsum(dd_june),cdd_july = cumsum(dd_july),cdd_aug = cumsum(dd_aug),cdd_early = cumsum(dd_early),cdd_mid = cumsum(dd_mid),cdd_ac = cumsum(dd_ac),cdd_mb = cumsum(dd_mb),cdd_ms = cumsum(dd_ms)

#load hopper data
setwd( paste(fdir, "grasshoppers/SexCombined/", sep="") )
hop= read.csv("HopperData_May2018.csv")

#======================================================
#CALCULATE GDD METRICS

#subset years to survey
clim1= clim[which(clim$Year %in% c(1958, 1959, 1960, 2006:2015))  ,]

#set temp outside summer to NA
clim1$Mean= (clim1$Min + clim1$Max)/2
clim1$Mean_summer= clim1$Mean
clim1[clim1$Julian<60 | clim1$Julian>243,"Mean_summer"]=NA

cdat=clim1

#metrics across years
clim1= ddply(clim1, c("Site", "Year"), summarise,
             Mean = mean(Mean_summer, na.rm=TRUE), Cdd_seas = max(cdd_sum, na.rm=FALSE), Cdd_seasfall = max(cdd_sumfall, na.rm=FALSE), Cdd_july = max(cdd_july, na.rm=TRUE) )
         #    Mean = mean(Mean_summer, na.rm=TRUE), Sd = sd(Mean_summer, na.rm=TRUE),Cdd_seas = max(cdd_sum, na.rm=TRUE),Cdd_june = max(cdd_june, na.rm=TRUE),Cdd_july = max(cdd_july, na.rm=TRUE),Cdd_aug = max(cdd_aug, na.rm=TRUE),Cdd_early = max(cdd_early, na.rm=TRUE),Cdd_mid = max(cdd_mid, na.rm=TRUE),Cdd_ac = max(cdd_ac, na.rm=TRUE),Cdd_mb = max(cdd_mb, na.rm=TRUE),Cdd_ms = max(cdd_ms, na.rm=TRUE), Sem = sd(Mean_summer, na.rm=TRUE)/sqrt(length(na.omit(Mean_summer))))

#PLOTS
#temp
ggplot(data=clim1, aes(x=Year, y = Mean, color=Site ))+geom_point()+geom_line() +theme_bw()
#cum dd
ggplot(data=clim1, aes(x=Year, y = Cdd_seas, color=Site ))+geom_point()+geom_line() +theme_bw()

#---------
# Jadult and GDDadult ~year by sites

#subset to dates with adults
hop1= hop[which(hop$in6>0),]

#subset to focal species
#specs= c("Aeropedellus clavatus","Chloealtis abdominalis","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes")

#hop1= hop1[which(hop1$species %in% specs ),]

#trim columns
hop1= hop1[,c("ordinal","species","in6","in5","in4","in3","in2","in1","total","site","period","year","sjy","cdd_sum")]

#-------------
## calculate median across individuals
hop1= hop1[order(hop1$ordinal),]
#cumulative sum of individuals within groups
hop1 = hop1 %>% group_by(species,site,year) %>% arrange(species,site,year,ordinal) %>% mutate(csind = cumsum(in6))
#number of median individual
hop3 = hop1 %>% group_by(species,site,year) %>% arrange(species,site,year,ordinal) %>% mutate(medind = max(csind)/2, qlowind=max(csind)*0.15, qupind=max(csind)*0.85)

#date of median individual
hop3$inddif= abs(hop3$medind-hop3$csind) #difference from median individual
hop3$inddif.qlow= abs(hop3$qlowind-hop3$csind) #difference from q20 individual
hop3$inddif.qup= abs(hop3$qupind-hop3$csind) #difference from q80 individual

hop4= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif),]))
hop4.qlow= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif.qlow),]))
hop4.qup= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif.qup),]))

#plot low and high percentiles
#combine
hop4$quantile=50
hop4.qup$quantile=85
hop4.qlow$quantile=15
hop5= rbind(hop4,hop4.qup, hop4.qlow)
hop5$quantile= as.factor(hop5$quantile)
hop5$yr_q= paste(hop5$year, hop5$quantile, sep="_")
hop5$year= as.factor(hop5$year)

hop4=hop5 #add quantiles 

#-----------------
# match phenology to temp and dd
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
hop4$siteyear= paste(hop4$site, hop4$year, sep="")

hop4$Tmean<-NA
match1= match(hop4$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop4$Tmean[matched]<- clim1$Mean[match1[matched]]  

hop4$cdd_seas<-NA
hop4$cdd_seas[matched]<- clim1$Cdd_seas[match1[matched]]  

#---------------
#add elevation
hop4$elevation= as.factor(elevs[match(hop4$site, sites)])

#---------------
#set up clim var

#GDDS
#add elevation
clim1$elevation= elevs[match(clim1$Site, sites)]
#add period
clim1$period="resurvey"
clim1$period[which(clim1$Year<2000)]<-"initial"

#order varaibles
#period
clim1$period= factor(clim1$period, levels=c("resurvey", "initial") )

#Order by average across sites with complete gdd data
clim.ave= subset(clim1, clim1$Site %in% c("B1","C1"))
clim.ave= aggregate(clim.ave, list(clim.ave$Year),FUN=mean )
clim1$Cdd_siteave= clim.ave$Cdd_seas[match(clim1$Year, clim.ave$Year)]
clim1$Cdd_july_siteave= clim.ave$Cdd_july[match(clim1$Year, clim.ave$Year)]

#------------------------------------------------
#CALCULATE DEVELOMENT INDEX

dat=hop[,c("ordinal","species","in6","in5","in4","in3","in2","in1","total","year","site","period","sjy","dd","dd_sum","cdd_sum","cdd_sumfall")]

#replace instar NAs with zero
dat[which(is.na(dat$in1)), "in1"]=0 
dat[which(is.na(dat$in2)), "in2"]=0
dat[which(is.na(dat$in3)), "in3"]=0
dat[which(is.na(dat$in4)), "in4"]=0
dat[which(is.na(dat$in5)), "in5"]=0
dat[which(is.na(dat$in6)), "in6"]=0

#Calculate development index
dat$DI=0
inds=which(dat$total>0)  
dat$DI[inds]= (dat$in1[inds] +dat$in2[inds]*2 +dat$in3[inds]*3 +dat$in4[inds]*4 +dat$in5[inds]*5 +dat$in6[inds]*6)/dat$total[inds]

#----------
# specs= c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus dawsoni", "Melanoplus sanguinipes", "Melanoplus boulderensis")
# #reduce to focal species
#dat= dat[dat$species %in% specs,] #[3:6]
#order
#dat$species= ordered(dat$species, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Chloealtis abdominalis", "Camnula pellucida","Melanoplus sanguinipes","Melanoplus dawsoni") )

#code period
dat$per=1
dat$per[dat$year>2000]=2

# #code early and late species
# dat$early_late=2
# dat$early_late[dat$species %in% specs[c(2,4)]]=1

#add elevation
dat$elev= as.factor(elevs[match(dat$site, sites)])

#add seasonal GDDs
dat$siteyear= paste(dat$site, dat$year, sep="")

dat$Tmean<-NA
match1= match(dat$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
dat$Tmean[matched]<- clim1$Mean[match1[matched]]  

dat$cdd_seas<-NA
dat$cdd_seas[matched]<- clim1$Cdd_seas[match1[matched]]  

dat$Cdd_siteave<-NA
dat$Cdd_siteave[matched]<- clim1$Cdd_siteave[match1[matched]]  

dat$Cdd_july_siteave<-NA
dat$Cdd_july_siteave[matched]<- clim1$Cdd_july_siteave[match1[matched]]  

#clean up
dat$year= as.factor(dat$year)

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#order varaibles
#period
dat$period= factor(dat$period, levels=c("resurvey", "initial") )
#elevation
dat$elev= factor(dat$elev, levels=c(3048,2591,2195,1752) )
#species
#dat$species= factor(dat$species, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Chloealtis abdominalis", "Camnula pellucida", "Melanoplus sanguinipes", "Melanoplus dawsoni") )

#------------------------------------------------
#ESTIMATE ADULTHOOD BASED ON DI

## failed using broom
#library(broom)
#broom::augment(x=fm1, newdata = Data, type.predict = "response")

#Indexed calculation
dat$spsiteyear= paste(dat$siteyear, dat$species, sep="")
combs= unique(dat$spsiteyear)

#days to predict over
doys= 150:265

#make matrix to store output
dout= data.frame(spsiteyear=combs, doy_adult= rep(NA, length(combs)),gdd_adult= rep(NA, length(combs)) ) 

for(k in 1:length(combs)){
  dats= subset(dat, dat$spsiteyear==combs[k])
  
  #require at least 5 data points
  if(nrow(dats)>=5) { 
    #doy
    doys= seq(min(dats$ordinal), max(dats$ordinal+7),5)
    
    spl<- smooth.spline(x=dats$ordinal, y=dats$DI)
    pred.spl<- predict(spl, doys)
    #extract point where almost all adults DI>5.5
    dout[k,2]= doys[which.max(pred.spl$y>5.5)]
    
    #gdd
    #restrict to observed gdds
    gdds= seq(min(dats$cdd_sumfall), max(dats$cdd_sumfall+50),10)
    
    spl<- smooth.spline(x=dats$cdd_sumfall, y=dats$DI)
    pred.spl<- predict(spl, gdds)
    #extract point where almost all adults DI>5.5
    dout[k,3]= gdds[which.max(pred.spl$y>5.5)]
    
  } #end check length
  
} #end combs

#add estimate back to df
dat$doy_adult= dout[match(dat$spsiteyear, dout$spsiteyear),"doy_adult"]
dat$gdd_adult= dout[match(dat$spsiteyear, dout$spsiteyear),"gdd_adult"]

#==============================================
#WRITE SPECIES AND CLIMATE DATA

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
write.csv(dat, "HopperClimateData.csv")
write.csv(clim1, "Clim1Data.csv")
#==============================================


