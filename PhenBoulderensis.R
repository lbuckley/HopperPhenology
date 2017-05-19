#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

#just Boulderensis
hop= hop[hop$species=="Melanoplus_boulderensis",]

#PLOTS
#Adult phenology
#Development index

#---------
#**************
#ADULT PHENOLOGY
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


#focus on sites B1 and C1 for now
hop3= hop2[which(hop2$site %in% c("B1", "C1", "CHA")) ,]

#drop due to missing climate data?
hop3= hop3[-which(hop3$siteyear %in% c("C12013","B12009")) ,] #

#add elevation
hop3$elevation= as.factor(elevs[match(hop3$site, sites)])

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#min ordinal date by dd, #plot regression
pdf("Phen_byGDD.pdf", height = 7, width = 10)
ggplot(data=hop3, aes(x=cdd, y = min, color=elevation))+geom_point(aes(shape=period, fill=period), size=3)+facet_wrap(~species, ncol=3) +theme_bw()+geom_smooth(method="lm")+ylim(125,250)+ylab("First appearance date")+xlab("Growing degree days")+ scale_color_manual(values=c("darkgreen","darkorange", "blue"))
dev.off()
