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

#fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#---------------------
#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd_sum = cumsum(dd_sum),cdd_june = cumsum(dd_june),cdd_july = cumsum(dd_july),cdd_aug = cumsum(dd_aug),cdd_early = cumsum(dd_early),cdd_mid = cumsum(dd_mid),cdd_ac = cumsum(dd_ac),cdd_mb = cumsum(dd_mb),cdd_ms = cumsum(dd_ms) ) 

#load hopper data
setwd( paste(fdir, "grasshoppers/SexCombined/", sep="") )
hop= read.csv("HopperData.csv")

#======================================================
#CALCULATE GDD METRICS

#subset years to survey
clim1= clim[which(clim$Year %in% c(1958, 1959, 1960, 2006:2015))  ,]

#set temp outside summer to NA
clim1$Mean= (clim1$Min + clim1$Max)/2
clim1$Mean_summer= clim1$Mean
clim1[clim1$Julian<60 | clim1$Julian>243,"Mean_summer"]=NA

#metrics across years
clim1= ddply(clim1, c("Site", "Year"), summarise,
             Mean = mean(Mean_summer, na.rm=FALSE), Cdd_seas = max(cdd_sum, na.rm=FALSE) )
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
specs= c("Aeropedellus clavatus","Chloealtis abdominalis","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes")

hop1= hop1[which(hop1$species %in% specs ),]

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

#==============================================
#==============================================
#FIGURE 1
#ELEVATION PLOTS

#GDDS
#add elevation
clim1$elevation= elevs[match(clim1$Site, sites)]
#add period
clim1$period="resurvey"
clim1$period[which(clim1$Year<2000)]<-"initial"

#Order by average across sites with complete gdd data
clim.ave= subset(clim1, clim1$Site %in% c("B1","C1"))
clim.ave= aggregate(clim.ave, list(clim.ave$Year),FUN=mean )
clim1$Cdd_siteave= clim.ave$Cdd_seas[match(clim1$Year, clim.ave$Year)]

plot_gdd_elev=ggplot(data=clim1, aes(x=elevation, y = Cdd_seas, group=Year, color=Cdd_siteave, linetype=period))+
  geom_line()+ 
  scale_colour_gradientn(name="Mean GDD", colours =matlab.like(10))+
  theme_bw()+ylab("season growing degree days (C)")+xlab("Elevation (m)")
#-------
#Plot ave phenology accross years

#ave across years
hop.ave= hop4[which(hop4$quantile==50),]
hop.ave$elevation= as.numeric(as.character(hop.ave$elevation))
hop.ave= as.data.frame( hop.ave[,c("elevation", "cdd_sum","site","year","ordinal","cdd_seas","species")] )
hop.ave= aggregate(hop.ave, list(hop.ave$elevation, hop.ave$species),FUN=mean )
hop.ave$species= hop.ave$Group.2

plot_gdds_elev=ggplot(data=hop.ave, aes(x=elevation, y = cdd_sum, group=species, color=species))+
  geom_line()+ 
  theme_bw()+ylab("cummulative growing degree days (C)")+xlab("elevation (m)")+theme(legend.position="none")

plot_doy_elev=ggplot(data=hop.ave, aes(x=elevation, y = ordinal, group=species, color=species))+
  geom_line()+ 
  theme_bw()+ylab("day of year")+xlab("elevation (m)")

#----------

#plot together
fig0= plot_grid(plot_gdd_elev, plot_gdds_elev, plot_doy_elev, labels = c('A', 'B','C'), rel_widths=c(1,0.7,1.3), nrow=1)

pdf("GDD_phen_byElev.pdf", height = 5, width = 12)
fig0
dev.off()
#===============================

#plot
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenology/figures/")

## FIGURE 2.
#Group by elevation rather than species
hop4q= subset(hop4, hop4$quantile==50)
#make elevation factor
hop4q$elevation= factor(hop4q$elevation, levels=c(3048,2591,2195,1752) )

#calculate significant regressions
#apply through combinations of species and elevations
hop4q$elevspec= paste(hop4q$elevation, hop4q$species, sep="")
elevspec= matrix(unique(hop4q$elevspec))
#extract p-values
p.doy= apply(elevspec,1, FUN=function(x) summary(lm(hop4q$ordinal[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)))] ~ hop4$cdd_seas[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(hop4q$cdd_sum[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)))] ~ hop4$cdd_seas[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
#combine
p.mat=as.data.frame(cbind(elevspec,p.doy,p.gdd))
#add columns for significance
p.mat$sig.doy<-"ns"
p.mat$sig.doy[which(p.doy<0.05)]="significant"
p.mat$sig.gdd="ns"
p.mat$sig.gdd[which(p.gdd<0.05)]="significant"

#add back to matrix
match1= match(hop4q$elevspec, elevspec)
hop4q$sig.doy= factor(p.mat[match1,"sig.doy"], levels=c("significant","ns"))
hop4q$sig.gdd= factor(p.mat[match1,"sig.gdd"], levels=c("significant","ns"))

#DOY
plot.phen=ggplot(data=hop4q, aes(x=cdd_seas, y = ordinal, color=species))+
  geom_point(aes(shape=period, fill=species, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, fill=NULL, stroke=1), size=3)+
  geom_smooth(method="lm",se=F, aes(linetype=sig.doy))+
  facet_wrap(~elevation, ncol=1, scales="free") +
  theme_bw()+ylab("day of year")+xlab("season growing degree days (C)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))+theme(legend.position="none")

#GDD
plot.phen.gdd=ggplot(data=hop4q, aes(x=cdd_seas, y = cdd_sum, color=species))+
  geom_point(aes(shape=period, fill=species, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, fill=NULL, stroke=1), size=3)+
  geom_smooth(method="lm",se=F, aes(linetype=sig.gdd))+
  facet_wrap(~elevation, ncol=1, scales="free") +
  theme_bw()+ylab("cummulative growing degree days")+xlab("season growing degree days (C)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))

#----
pdf("plot_phen.pdf", height = 12, width = 10)
plot_grid(plot.phen, plot.phen.gdd, nrow=1, rel_widths=c(1,1.5) )
dev.off()

#====================================
## FIGURE 3.
# DEVELOPMENT INDEX

dat=hop[,c("ordinal","species","in6","in5","in4","in3","in2","in1","total","year","site","period","sjy","cdd_sum")]

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

#ANALYSIS
specs= c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus dawsoni", "Melanoplus sanguinipes", "Melanoplus boulderensis")

#reduce to focal species
dat= dat[dat$species %in% specs[2:5],]

#code period
dat$per=1
dat$per[dat$year>2000]=2

#code early and late species
dat$early_late=2
dat$early_late[dat$species %in% specs[c(2,4)]]=1

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

#------------------------------------------------------------
#PLOTS
dat$species= as.factor(dat$species)
dat$year= as.factor(dat$year)

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#--------------------
#DEVELOPMENTAL INDEX
#Plot DI by ordinal date
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=cdd_seas, group=siteyear, linetype=period))+facet_grid(species~elev) +
  theme_bw()+geom_smooth(se=FALSE)+
  scale_colour_gradientn(colours =matlab.like(10))
#***** 

#Plot DI by GDD
dat$cdd= dat$cdd_sum
di.plot= ggplot(data=dat, aes(x=cdd_sum, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+geom_line()+xlim(0,600)
#note xlim restricted

#==============================================================================
##Bin by GDD
gdds= seq( min(dat$cdd_sum, na.rm=TRUE), max(dat$cdd_sum, na.rm=TRUE), length.out=20) #CHANGE NUMBER GDD BINS
dat$gdd.bin= cut(dat$cdd_sum, breaks = gdds, labels=FALSE)

dat.gddbin = dat %>% group_by(species,site,year,gdd.bin) %>% summarise_each(funs(mean))

##Bin by date
dates= seq( min(dat$ordinal, na.rm=TRUE)-1, max(dat$ordinal, na.rm=TRUE)+1, 20)
dat$date.bin= cut(dat$ordinal, breaks = dates, labels=FALSE, include.highest=TRUE)

dat.datebin = dat %>% group_by(species,site,year,date.bin) %>% summarise_each(funs(mean))

#----------------
#group by seasonal GDD
#drop sites without data
clim.seas2= clim1[which(clim1$Site %in% c("B1","C1")),]
clim.seas2= clim.seas2 %>% group_by(Year) %>% summarise(Cdd_seas= mean(Cdd_seas) )
clim.seas2= clim.seas2[order(clim.seas2$Cdd_seas),]

#------
#CODE BY GDD
#initial 1958,1959,1960
#cold 2009,2010,2014
#med 2008,2011,2013,2015
#warm 2006, 2007,2012

dat$per=NA
dat$per[which(dat$year %in% c(1958,1959,1960) )]="initial"
dat$per[which(dat$year %in% c(2009,2010,2014) )]="cold"
dat$per[which(dat$year %in% c(2008,2011,2013,2015) )]="med"
dat$per[which(dat$year %in% c(2006,2007,2012) )]="warm"

#--------------------------
#PLOT DI by grouped years

#average across years
dat$gdd.binned= gdds[dat$gdd.bin]
dat2= dat %>% group_by(species,elev, site, per,gdd.binned) %>%  mutate(DI=mean(DI) )

#order varaibles
#period
dat$per= factor(dat$per, levels=c("initial","cold","med","warm") )
#elevation
dat$elev= factor(dat$elev, levels=c(3048,2591,2195,1752) )
#species
dat$species= factor(dat$species, levels=c("Melanoplus boulderensis", "Camnula pellucida", "Melanoplus sanguinipes", "Melanoplus dawsoni") )

#DEVELOPMENTAL INDEX
#Plot DI by ordinal date
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=per))+facet_grid(elev~species) +
  theme_bw()+geom_smooth(se=F) +
  scale_color_manual(values=c("black", "blue", "green","red") ) +
  theme_bw()+ylab("development index")+xlab("day of year")+labs(color="season")
#scale_linetype_manual(values=c("dashed", "solid", "solid","solid") ) +
#+ geom_point(aes(shape=period), size=2)
## USE FIG? *************

ggplot(data=dat, aes(x=ordinal, y = DI, color=year))+facet_grid(elev~species) +
  theme_bw()+geom_smooth(se=F) +
  scale_color_manual(values=c("black", "blue", "green","red") ) +
  theme_bw()+ylab("development index")+xlab("day of year")+labs(color="season")


