library(reshape2)

#LOAD data
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/grasshoppers/ForLTER/")

#-----------------
#juveniles
b1.j= read.csv("HopperData_Sep2019_NymphalB1_juvenile.csv")
c1.j= read.csv("HopperData_Sep2019_NymphalC1_juvenile.csv")
chaut.j= read.csv("HopperData_Sep2019_NymphalChaut_juvenile.csv")
#add sites
b1.j$site="B1"
c1.j$site="C1"
chaut.j$site="CHA"
#bind
nj= rbind(b1.j, c1.j, chaut.j)
#add date
date= as.Date(nj$date, format="%m/%d/%y")
nj$ordinal= format(date,"%j")
nj$ordinal= as.integer(nj$ordinal)
nj$year= format(date,"%Y")
nj$sp.site.yr.ord= paste(nj$Species, nj$site, nj$year, nj$ordinal, sep="_")

#-----------------
#adults
b1.a= t(read.csv("HopperData_Sep2019_NymphalB1_adult.csv", header=FALSE) )
c1.a= t(read.csv("HopperData_Sep2019_NymphalC1_adult.csv", header=FALSE) )
chaut.a= t(read.csv("HopperData_Sep2019_NymphalChaut_adult.csv", header=FALSE) )
#fix format
colnames(b1.a)<- b1.a[1,]
colnames(c1.a)<- c1.a[1,]
colnames(chaut.a)<- chaut.a[1,]
b1.a=b1.a[-1,]
c1.a=c1.a[-1,]
chaut.a=chaut.a[-1,]
b1.a= as.data.frame(b1.a)
c1.a= as.data.frame(c1.a)
chaut.a= as.data.frame(chaut.a)
#to long format
b1.a.m= melt(b1.a, id = c("DATE","Ordinal Date:"))
c1.a.m= melt(c1.a, id = c("DATE","OrdinalDate"))
chaut.a.m= melt(chaut.a, id = c("Date","Ordinal Date:"))
colnames(b1.a.m)<- colnames(c1.a.m)
colnames(chaut.a.m)<- colnames(c1.a.m)
#add sites
b1.a.m$site="B1"
c1.a.m$site="C1"
chaut.a.m$site="CHA"
#bind
nad= rbind(b1.a.m, c1.a.m, chaut.a.m)

#add date
date= as.Date(nad$DATE, format="%m/%d/%y")
nad$ordinal= format(date,"%j")
nad$ordinal= as.integer(nad$ordinal)
nad$year= format(date,"%Y")
nad$sp.site.yr.ord= paste(nad$variable, nad$site, nad$year, nad$ordinal, sep="_")

#match format
colnames(nj)[6]<-"X6.Adult"
nad$GDDs<-NA
nad1= cbind(nad[,c("year","DATE","ordinal","GDDs","variable","value")], matrix(NA, nrow(nad),5), nad[,c("site","sp.site.yr.ord")] )   
colnames(nad1)<- colnames(nj)

#--------------
#match
nad1$X6.Adult= as.numeric(as.character(nad1$X6.Adult))
nj$X6.Adult= NA
nj$X6.Adult= as.numeric(nj$X6.Adult)

match1=match(nj$sp.site.yr.ord, nad1$sp.site.yr.ord)
matched= which(!is.na(match1))

nj$X6.Adult[matched]= nad1$X6.Adult[match1[matched]] 

#add unmatched
match2=match(nad1$sp.site.yr.ord, nj$sp.site.yr.ord)
unmatched= which(is.na(match2))

nall= rbind(nj, nad1[unmatched,])

date= as.Date(nall$date, format="%m/%d/%y")
nall$month= format(date,"%m")
nall$month= as.integer(nall$month)
nall$day= format(date,"%d")
nall$day= as.integer(nall$day)
nall$period<- "resurvey"

#-----------------
#add historic
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/grasshoppers/ForLTER/")
hj= read.csv("HopperData_NymphalHist.csv")

hj$DateCollected= format(as.Date(hj$DateCollected, "%m/%d/%y"), "19%y-%m-%d")

date= as.Date(hj$DateCollected, format="%Y-%m-%d")
hj$ordinal= format(date,"%j")
hj$ordinal= as.integer(hj$ordinal)
hj$year= format(date,"%Y")
hj$month= format(date,"%m")
hj$month= as.integer(hj$month)
hj$day= format(date,"%d")
hj$day= as.integer(hj$day)
hj$period<- "initial"

hj$GDDs= NA
#fix site
hj$Site= as.character(hj$Site)
hj$Site[which(hj$Site=="Chaut")]<-"CHA"
hj$sp.site.yr.ord= paste(hj$GenSpec, hj$Site, hj$year, hj$ordinal, sep="_")

hj.bind= hj[,c("year", "DateCollected", "ordinal", "GDDs", "GenSpec", "SumAdult", "SumN5", "SumN4", "SumN3", "SumN2", "SumN1", "Site", "sp.site.yr.ord","month","day","period")]    
names(hj.bind)<-names(nall) 
#combine
nall= rbind(nall, hj.bind) 
  
#--------------
#WRITE OUT
write.csv(nall, "NymphalData.csv")

#=========================================
#ADD CLIMATE DATA

fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled_May2018.csv")

#load hopper data
setwd( paste(fdir, "grasshoppers/SexCombined/", sep="") )
hop= read.csv("HopperData_May2018.csv")

#------------------------------------
#Add GDD data

#Change to NOAA for matching
nall[which(nall$site=="CHA"),"site"]<-"NOAA"
#columns for matching
nall$sjy= paste(nall$site, nall$ordinal, nall$year, sep="_")

#match
match1= match(nall$sjy,clim$sjy)
matched= which(!is.na(match1))

nall$dd=NA;nall$dd_sum=NA;nall$cdd_sum=NA; nall$cdd_sumfall=NA;
nall[matched,c("dd","dd_sum","cdd_sum","cdd_sumfall")]= clim[match1[matched], c("dd","dd_sum","cdd_sum","cdd_sumfall")]

#--------
#combine with other hopper data
hop=hop[,-1]

nall$total= rowSums(nall[,6:11], na.rm=TRUE)
nall.bind= nall[,c("date", "ordinal", "GDDs", "Species", "X6.Adult", "X5", "X4", "X3", "X2", "X1","total","month", "day","year","site","period", "sjy", "dd", "dd_sum", "cdd_sum", "cdd_sumfall")]
names(nall.bind)<- names(hop)
hop= rbind(hop, nall.bind)
#fix NOAA sites
hop$site[which(hop$site=="NOAA")]<-"CHA"

#Write hopper data out
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/grasshoppers/ForLTER/final/")
write.csv(hop, "HopperData_Sept2019.csv")

#====================================================

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

setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/grasshoppers/ForLTER/final/")
hop= read.csv("HopperData_Sept2019.csv")

#======================================================
#READ DATA

#fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled_May2018.csv")

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
dat$species= factor(dat$species) #, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Chloealtis abdominalis", "Camnula pellucida", "Melanoplus sanguinipes", "Melanoplus dawsoni") )

#add wind length
dat$species.lab= paste(dat$species,"(SW)", sep=" ")
#fix long wing
dat$species.lab= replace(dat$species.lab, which(dat$species.lab=="Camnula pellucida (SW)"), "Camnula pellucida (LW)")
dat$species.lab= replace(dat$species.lab, which(dat$species.lab=="Melanoplus sanguinipes (SW)"), "Melanoplus sanguinipes (LW)")
dat$species.lab= factor(dat$species.lab, levels=c("Aeropedellus clavatus (SW)","Melanoplus boulderensis (SW)","Chloealtis abdominalis (SW)", "Camnula pellucida (LW)", "Melanoplus sanguinipes (LW)", "Melanoplus dawsoni (SW)") )

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

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
write.csv(dat, "HopperData_Sept2019_forPhenOverlap.csv")

#======================================================
