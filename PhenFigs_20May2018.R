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
             Mean = mean(Mean_summer, na.rm=TRUE), Sd = sd(Mean_summer, na.rm=TRUE),Cdd_sum = max(cdd_sum, na.rm=TRUE),Cdd_june = max(cdd_june, na.rm=TRUE),Cdd_july = max(cdd_july, na.rm=TRUE),Cdd_aug = max(cdd_aug, na.rm=TRUE),Cdd_early = max(cdd_early, na.rm=TRUE),Cdd_mid = max(cdd_mid, na.rm=TRUE),Cdd_ac = max(cdd_ac, na.rm=TRUE),Cdd_mb = max(cdd_mb, na.rm=TRUE),Cdd_ms = max(cdd_ms, na.rm=TRUE), Sem = sd(Mean_summer, na.rm=TRUE)/sqrt(length(na.omit(Mean_summer))))

#PLOTS
#temp
ggplot(data=clim1, aes(x=Year, y = Mean, color=Site ))+geom_point()+geom_line() +theme_bw()
#cum dd
ggplot(data=clim1, aes(x=Year, y = Cdd_sum, color=Site ))+geom_point()+geom_line() +theme_bw()

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

#order species by seasonal timing
#ave phen
#hop.el = hop1 %>% group_by(species,site) %>% arrange(species,site) %>% mutate(phen = mean(ordinal))
hop.agg= aggregate(hop1, list(hop1$species),FUN=mean) #fix for hop1$site,
hop.agg= hop.agg[order(hop.agg$ordinal),]

#make species factor for plotting
hop1$species= factor(hop1$species, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Camnula pellucida","Melanoplus sanguinipes", "Melanoplus dawsoni", "Chloealtis abdominalis"))
hop2$species= factor(hop2$species, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Camnula pellucida","Melanoplus sanguinipes", "Melanoplus dawsoni", "Chloealtis abdominalis"))

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
hop2$siteyear= paste(hop2$site, hop2$year, sep="")
hop4$siteyear= paste(hop4$site, hop4$year, sep="")

hop2$Tmean=NA
match1= match(hop2$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop2$Tmean[matched]<- clim1$Mean[match1[matched]]  

hop2$cdd_yr=NA
hop2$cdd_yr[matched]<- clim1$Cdd_sum[match1[matched]]  

#----------------------------------------------
hop4$Tmean=NA
match1= match(hop4$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop4$Tmean[matched]<- clim1$Mean[match1[matched]]  

hop4$cdd_yr=NA
hop4$cdd_yr[matched]<- clim1$Cdd_sum[match1[matched]]  

#----------------------------------------------

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

##Specify cdd metric
hop2$cdd=NA; hop4$cdd=NA
hop2$cdd=hop2$cdd_sum; hop4$cdd=hop4$cdd_sum

#min ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = min, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = mean, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#median individual ordinal by dd
ggplot(data=hop4, aes(x=cdd, y = ordinal, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
  
#drop sites with missing gdd data
hop2= hop2[which(hop2$cdd_sum>0),]
hop4= hop4[which(hop4$cdd_sum>0),]

#add elevation
hop2$elevation= as.factor(elevs[match(hop2$site, sites)])
hop4$elevation= as.factor(elevs[match(hop4$site, sites)])


#===============================
#FIGURE 1
#ELEVATION PLOTS

#GDDS
#add elevation
clim1$elevation= elevs[match(clim1$Site, sites)]
#fix NAs
clim1$Cdd_sum[which(is.infinite(clim1$Cdd_sum))]<- NA
#add period
clim1$period="resurvey"
clim1$period[which(clim1$Year<2000)]<-"initial"

#Order by average across sites with complete gdd data
clim.ave= subset(clim1, clim1$Site %in% c("B1","C1"))
clim.ave= aggregate(clim.ave, list(clim.ave$Year),FUN=mean )
clim1$Cdd_siteave= clim.ave$Cdd_sum[match(clim1$Year, clim.ave$Year)]

plot_gdd_elev=ggplot(data=clim1, aes(x=elevation, y = Cdd_sum, group=Year, color=Cdd_siteave, linetype=period))+
  geom_line()+ 
  scale_colour_gradientn(name="Mean GDD", colours =matlab.like(10))+
  theme_bw()+ylab("season growing degree days (C)")+xlab("Elevation (m)")
#-------
#Plot ave phenology accross years

#ave across years
hop.ave= hop4[which(hop4$quantile==50),]
hop.ave$elevation= as.numeric(as.character(hop.ave$elevation))
hop.ave= as.data.frame( hop.ave[,c("elevation", "cdd_sum","site","year","ordinal","cdd","species")] )
hop.ave= aggregate(hop.ave, list(hop.ave$elevation, hop.ave$species),FUN=mean )
hop.ave$species= hop.ave$Group.2

plot_gdds_elev=ggplot(data=hop.ave, aes(x=elevation, y = cdd, group=species, color=species))+
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
p.doy= apply(elevspec,1, FUN=function(x) summary(lm(hop4q$ordinal[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)))] ~ hop4$cdd_yr[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(hop4q$cdd_yr[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)))] ~ hop4$cdd_yr[which(hop4q$elevation==substr(x,1,4)&hop4q$species==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
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
plot.phen=ggplot(data=hop4q, aes(x=cdd_yr, y = ordinal, color=species))+
  geom_point(aes(shape=period, fill=species, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, fill=NULL, stroke=1), size=3)+
  geom_smooth(method="lm",se=F, aes(linetype=sig.doy))+
  facet_wrap(~elevation, ncol=1, scales="free") +
  theme_bw()+ylab("day of year")+xlab("season growing degree days (C)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))+theme(legend.position="none")

#GDD
plot.phen.gdd=ggplot(data=hop4q, aes(x=cdd_yr, y = cdd, color=species))+
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

#===============
