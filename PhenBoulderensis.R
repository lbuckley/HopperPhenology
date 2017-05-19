#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

sites= c("Redfox", "A1", "B1", "C1", "D1")  
elevs= c(1574, 2195, 2591, 3048, 3739)

#just Boulderensis
hop1= hop[hop$species=="Melanoplus boulderensis",]

#PLOTS
#Adult phenology
#Development index

#---------
#ADULT PHENOLOGY
# Jadult and GDDadult ~year by sites

#subset to dates with adults
hop1= hop1[which(hop1$in6>0),]

#metrics across years
hop2= ddply(hop1, c("site", "year","period","species"), summarise,
            min = min(ordinal, na.rm=TRUE), mean = mean(ordinal, na.rm=TRUE) )

#Change years for plotting ### FIX
hop2[which(hop2$year==1958),"year"]= 1958+40
hop2[which(hop2$year==1959),"year"]= 1959+40
hop2[which(hop2$year==1960),"year"]= 1960+40

#min ordinal date
ggplot(data=hop2, aes(x=year, y = min, color=site, shape=period))+geom_point()+geom_line() +theme_bw()
#mean ordinal date
ggplot(data=hop2, aes(x=year, y = mean, color=site, shape=period))+geom_point()+geom_line()+theme_bw()

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
ggplot(data=hop4, aes(x=year, y = ordinal, color=site ))+geom_point()+geom_line() +theme_bw()

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

#Calculate cdd in month before mean ordinal date
hop2$cdd_ss=NA
hop2$cdd_ss[matched]<- clim1$Cdd_mb[match1[matched]] 
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

#Calculate cdd in month before mean ordinal date
hop4$cdd_ss=NA
hop4$cdd_ss[matched]<- clim1$Cdd_mb[match1[matched]]

#----------------------------------------------

#TEMP PLOTS
#min ordinal date by temp
hop2= hop2[order(hop2$Tmean),]
ggplot(data=hop2, aes(x=Tmean, y = min, color=site,shape=period))+geom_point()+geom_line() +theme_bw()
#mean ordinal date by temp
ggplot(data=hop2, aes(x=Tmean, y = mean, color=site,shape=period))+geom_point()+geom_line()+theme_bw()
#median individual ordinal by temp
ggplot(data=hop4, aes(x=Tmean, y = ordinal, color=site,shape=period))+geom_point()+geom_line()+theme_bw()

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
ggplot(data=hop2, aes(x=cdd, y = min, color=site, shape=period))+geom_point()+geom_line()+theme_bw()
#mean ordinal date by dd
ggplot(data=hop2, aes(x=cdd, y = mean, color=site, shape=period))+geom_point()+geom_line()+theme_bw()
#median individual ordinal by dd
ggplot(data=hop4, aes(x=cdd, y = ordinal, color=site, shape=period))+geom_point()+geom_line()+theme_bw()

hop3=hop2
#drop due to missing climate data?
hop3= hop3[-which(hop3$siteyear %in% c("C12013","B12009")) ,] #

#add elevation
hop3$elevation= as.factor(elevs[match(hop3$site, sites)])

#PLOTS
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#min ordinal date by dd, #plot regression
pdf("Phen_MBould.pdf", height = 7, width = 10)
ggplot(data=hop3, aes(x=cdd, y = min, color=elevation))+geom_point(aes(shape=period, fill=period), size=3)+theme_bw()+geom_smooth(method="lm")+ylim(125,250)+ylab("First appearance date")+xlab("Growing degree days")+ scale_color_manual(values=c("darkgreen","darkorange", "blue"))
dev.off()
