#load libraries
library(ggplot2)
library(plyr) 
library(dplyr)

sites= c("CHA", "A1", "B1", "C1", "D1")  #Redfox: 1574
elevs= c(1752, 2195, 2591, 3048, 3739)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

#======================================================
#READ DATA

fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#---------------------
#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd_sum = cumsum(dd_sum),cdd_june = cumsum(dd_june),cdd_july = cumsum(dd_july),cdd_aug = cumsum(dd_aug),cdd_early = cumsum(dd_early),cdd_mid = cumsum(dd_mid),cdd_ac = cumsum(dd_ac),cdd_mb = cumsum(dd_mb),cdd_ms = cumsum(dd_ms) ) 


#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )
hop= read.csv("HopperData.csv")

#======================================================
#ANALYSIS

#temperature plots
# temp, var, and GDD ~year by sites

#subset years to survey
clim1= clim[which(clim$Year %in% c(1958, 1959, 1960, 2006:2015))  ,]

##Change years for plotting ### FIX
#clim1[which(clim1$Year==1958),"Year"]= 1958+40
#clim1[which(clim1$Year==1959),"Year"]= 1959+40
#clim1[which(clim1$Year==1960),"Year"]= 1960+40

#set temp outside summer to NA
clim1$Mean= (clim1$Min + clim1$Max)/2
clim1$Mean_summer= clim1$Mean
clim1[clim1$Julian<60 | clim1$Julian>243,"Mean_summer"]=NA

#metrics across years
clim1= ddply(clim1, c("Site", "Year"), summarise,
             Mean = mean(Mean_summer, na.rm=TRUE), Sd = sd(Mean_summer, na.rm=TRUE),Cdd_sum = max(cdd_sum, na.rm=TRUE),Cdd_june = max(cdd_june, na.rm=TRUE),Cdd_july = max(cdd_july, na.rm=TRUE),Cdd_aug = max(cdd_aug, na.rm=TRUE),Cdd_early = max(cdd_early, na.rm=TRUE),Cdd_mid = max(cdd_mid, na.rm=TRUE),Cdd_ac = max(cdd_ac, na.rm=TRUE),Cdd_mb = max(cdd_mb, na.rm=TRUE),Cdd_ms = max(cdd_ms, na.rm=TRUE), Sem = sd(Mean_summer, na.rm=TRUE)/sqrt(length(na.omit(Mean_summer))))
#better NA handling?

#temp
ggplot(data=clim1, aes(x=Year, y = Mean, color=Site ))+geom_point()+geom_line() +theme_bw()
#sd
ggplot(data=clim1, aes(x=Year, y = Sd, color=Site ))+geom_point()+geom_line() +theme_bw()
#cum dd
ggplot(data=clim1, aes(x=Year, y = Cdd_sum, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_june, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_july, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_aug, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_early, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_mid, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_ac, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_mb, color=Site ))+geom_point()+geom_line() +theme_bw()
ggplot(data=clim1, aes(x=Year, y = Cdd_ms, color=Site ))+geom_point()+geom_line() +theme_bw()

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

##Change years for plotting ### FIX
#hop2[which(hop2$year==1958),"year"]= 1958+40
#hop2[which(hop2$year==1959),"year"]= 1959+40
#hop2[which(hop2$year==1960),"year"]= 1960+40

#order species by seasonal timing
#ave phen
#hop.el = hop1 %>% group_by(species,site) %>% arrange(species,site) %>% mutate(phen = mean(ordinal))
hop.agg= aggregate(hop1, list(hop1$species),FUN=mean) #fix for hop1$site,
hop.agg= hop.agg[order(hop.agg$ordinal),]

#make species factor for plotting
#hop1$species= factor(hop1$species, levels=hop.agg$Group.1)
#hop2$species= factor(hop2$species, levels=hop.agg$Group.1)
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
## refine dd metrics? restrict to summer?

clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
hop2$siteyear= paste(hop2$site, hop2$year, sep="")
hop4$siteyear= paste(hop4$site, hop4$year, sep="")

hop2$Tmean=NA
match1= match(hop2$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop2$Tmean[matched]<- clim1$Mean[match1[matched]]  

hop2$cdd_sum=NA
hop2$cdd_sum[matched]<- clim1$Cdd_sum[match1[matched]]  
hop2$cdd_june=NA
hop2$cdd_june[matched]<- clim1$Cdd_june[match1[matched]] 
hop2$cdd_july=NA
hop2$cdd_july[matched]<- clim1$Cdd_july[match1[matched]] 
hop2$cdd_aug=NA
hop2$cdd_aug[matched]<- clim1$Cdd_aug[match1[matched]] 
hop2$cdd_early=NA
hop2$cdd_early[matched]<- clim1$Cdd_early[match1[matched]] 
hop2$cdd_mid=NA
hop2$cdd_mid[matched]<- clim1$Cdd_mid[match1[matched]] 

#--------------------------
#Calculate cdd in month before mean ordinal date

hop2$cdd_ss=NA
hop2$cdd_ss[matched]<- clim1$Cdd_ms[match1[matched]] 

#clav
inds= which(hop2$species[matched]=="Aeropedellus clavatus")
hop2$cdd_ss[matched[inds]]<- clim1$Cdd_ac[match1[matched[inds]]] 

#bould
inds= which(hop2$species[matched]=="Melanoplus boulderensis")
hop2$cdd_ss[matched[inds]]<- clim1$Cdd_mb[match1[matched[inds]]] 

#----------------------------------------------
hop4$Tmean=NA
match1= match(hop4$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
hop4$Tmean[matched]<- clim1$Mean[match1[matched]]  

hop4$cdd_sum=NA
hop4$cdd_sum[matched]<- clim1$Cdd_sum[match1[matched]]  
hop4$cdd_june=NA
hop4$cdd_june[matched]<- clim1$Cdd_june[match1[matched]] 
hop4$cdd_july=NA
hop4$cdd_july[matched]<- clim1$Cdd_july[match1[matched]] 
hop4$cdd_aug=NA
hop4$cdd_aug[matched]<- clim1$Cdd_aug[match1[matched]] 
hop4$cdd_early=NA
hop4$cdd_early[matched]<- clim1$Cdd_early[match1[matched]] 
hop4$cdd_mid=NA
hop4$cdd_mid[matched]<- clim1$Cdd_mid[match1[matched]] 

#--------------------------
#Calculate cdd in month before mean ordinal date

hop4$cdd_ss=NA
hop4$cdd_ss[matched]<- clim1$Cdd_ms[match1[matched]] 

#clav
inds= which(hop4$species[matched]=="Aeropedellus clavatus")
hop4$cdd_ss[matched[inds]]<- clim1$Cdd_ac[match1[matched[inds]]] 

#bould
inds= which(hop4$species[matched]=="Melanoplus boulderensis")
hop4$cdd_ss[matched[inds]]<- clim1$Cdd_mb[match1[matched[inds]]] 

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
hop2$cdd=hop2$cdd_june; hop4$cdd=hop4$cdd_june
hop2$cdd=hop2$cdd_july; hop4$cdd=hop4$cdd_july
hop2$cdd=hop2$cdd_aug; hop4$cdd=hop4$cdd_aug
hop2$cdd=hop2$cdd_early; hop4$cdd=hop4$cdd_early
hop2$cdd=hop2$cdd_mid; hop4$cdd=hop4$cdd_mid
hop2$cdd=hop2$cdd_ss; hop4$cdd=hop4$cdd_ss

#min ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = min, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = mean, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#median individual ordinal by dd
ggplot(data=hop4, aes(x=cdd, y = ordinal, color=site, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()

plot.qs=ggplot(data=hop4, aes(x=cdd, y = ordinal, color=site, linetype=quantile))+geom_point(aes(shape=period), size=3)+geom_smooth(method="lm",se=F)+facet_wrap(~species, ncol=3) +theme_bw()

#drop sites with missing gdd data
hop2= hop2[which(hop2$cdd_sum>0),]

#focus on sites B1 and C1 for now
hop3= hop2[which(hop2$site %in% c("A1","B1", "C1", "CHA")) ,]

#add elevation
hop3$elevation= as.factor(elevs[match(hop3$site, sites)])

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenology\\figures\\")

#quantiles
pdf("PhenQuantiles_byGDD.pdf", height = 12, width = 12)
plot.qs
dev.off()

#min ordinal date by dd, #plot regression
pdf("Phen_byGDD.pdf", height = 7, width = 10)
ggplot(data=hop3, aes(x=cdd, y = min, color=elevation))+geom_point(aes(shape=period, fill=period), size=3)+facet_wrap(~species, ncol=3) +theme_bw()+geom_smooth(method="lm")+ylim(125,250)+ylab("First appearance date")+xlab("Growing degree days")+ scale_color_manual(values=c("darkgreen","darkorange", "blue","purple"))
dev.off()

#====================================
#ADDITIONAL GDD PLOTS

#subset to years
clim2= clim[which(clim$Year %in% c(1958:1960,2006:2015)),]
#subset sites
clim2= clim2[which(clim2$Site %in% c("A1","B1","C1","CHA")),]

#add time period
clim2$per=NA
clim2$per[which(clim2$Year %in% c(1958,1959,1960) )]="initial"
clim2$per[which(clim2$Year %in% c(2009,2010,2014) )]="cold"
clim2$per[which(clim2$Year %in% c(2008,2011,2013,2015) )]="med"
clim2$per[which(clim2$Year %in% c(2006,2007,2012) )]="warm"

clim2$cdd= clim2$cdd_sum

#plot GDD accumulation over time
gdd.plot= ggplot(data=clim2, aes(x=Julian, y = cdd, color=as.factor(Year),linetype=per))+geom_line()+facet_grid(.~Site, scales="free_y") +theme_bw()

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#min ordinal date by dd, #plot regression
pdf("GDD_plot.pdf", height = 7, width = 10)
gdd.plot
dev.off()

#------------------------------
#plot initial, recent GDDs

hop2$min.cdd=NA
hop2$mean.cdd=NA

clim$cdd= clim$cdd_sum

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
ggplot(data=hop2, aes(x=year, y = min.cdd, color=species))+geom_point()+geom_line()+facet_grid(.~site) +theme_bw()+ theme(legend.position = "bottom")
#gdd at peak
ggplot(data=hop2, aes(x=year, y = mean.cdd, color=species))+geom_point()+geom_line()+facet_grid(.~site) +theme_bw()+ theme(legend.position = "bottom")
