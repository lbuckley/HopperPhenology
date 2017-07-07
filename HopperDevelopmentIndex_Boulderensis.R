#load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(gridExtra)

library(lme4)
library(nlme)
library(grDevices)
library(akima) #for interpolation
library(grid)
library(gridExtra)

#Elev for A1  B1  C1  CHA
sites= c("A1", "B1", "C1", "CHA")
elevs= c(2195, 2591, 3048, 1752)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

count= function(x){
  length( !is.na(x))
}

#--------------------------------------
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )
hop= read.csv("HopperData.csv")
#fix GDD column
hop$GDDs= hop$cdd_sum
#-------------------------
#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd_sum = cumsum(dd_sum),cdd_june = cumsum(dd_june),cdd_july = cumsum(dd_july),cdd_aug = cumsum(dd_aug),cdd_early = cumsum(dd_early),cdd_mid = cumsum(dd_mid),cdd_ac = cumsum(dd_ac),cdd_mb = cumsum(dd_mb),cdd_ms = cumsum(dd_ms) )

#==========================================================

dat=hop

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

#In relation to OD
mod1= lme(DI ~ poly(ordinal,3)+elev*per* species, random=~1|year, data=dat)
anova(mod1)

#GDD TEXT MODEL
dat1=na.omit(dat)
#nlme model
mod1= lme(DI ~ poly(GDDs,3)+elev*per+species, random=~1|year, data=dat1)
mod1= lme(DI ~ poly(GDDs,3)+elev*per*species, random=~1|year, data=dat1)
anova(mod1)

dat1= subset(dat1, dat1$species==specs[[2]])
mod1= lme(DI ~ poly(GDDs,3)+elev*per, random=~1|year, data=dat1)

anova(mod1,type="marginal")

#------------------------------------------------------------
#PLOTS

dat$species= as.factor(dat$species)
dat$year= as.factor(dat$year)

#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\DevelopInd\\")
#file<-paste("DIbinned_ByGDD.pdf" ,sep="", collapse=NULL)
#pdf(file,height = 8, width = 11)

#par(mfrow=c(2,2), cex=1.2, mar=c(1.5, 2, 2, 0), mgp=c(2, 1, 0), oma=c(2,2,0,0), lwd=2)

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#--------------------
#DEVELOPMENTAL INDEX
#Plot DI by ordinal date
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+geom_line()

#Plot DI by GDD
dat$cdd= dat$cdd_sum
di.plot= ggplot(data=dat, aes(x=cdd, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+geom_line()+xlim(0,600)
#note xlim restricted

#==============================================================================
##Bin by GDD
gdds= seq( min(dat$GDDs, na.rm=TRUE), max(dat$GDDs, na.rm=TRUE), length.out=30)
dat$gdd.bin= cut(dat$GDD, breaks = gdds, labels=FALSE)

dat.gddbin = dat %>% group_by(species,site,year,gdd.bin) %>% summarise_each(funs(mean))

##Bin by date
dates= seq( min(dat$ordinal, na.rm=TRUE)-1, max(dat$ordinal, na.rm=TRUE)+1, 20)
dat$date.bin= cut(dat$ordinal, breaks = dates, labels=FALSE, include.highest=TRUE)

dat.datebin = dat %>% group_by(species,site,year,date.bin) %>% summarise_each(funs(mean))

#dat$timeper= cut(dat$year, breaks = c(1957,1961,2011), labels=FALSE)

#=======================================
#Just boulderensis

dat.t3= subset(dat.t1, dat.t1$species=="Melanoplus boulderensis")
dat.t3$GDDs_binned= gdds[dat.t3$gdd.bin]

#--------------

g1= ggplot(data=dat.t3) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(site~year, scales = "free_x")+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,300)

#--------------------
#split by site
dat.b1= subset(dat.t3, dat.t3$site %in% c("B1") )

g.b1= ggplot(data=dat.b1) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(.~year)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,300)

#by GDD
ggplot(data=dat.b1) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(.~dd.seas)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,300)

#--------------------
dat.c1= subset(dat.t3, dat.t3$site %in% c("C1") )

g.c1= ggplot(data=dat.c1) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(.~year)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,150)

#by GDD
ggplot(data=dat.c1) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(.~dd.seas)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,150)

#------------------
# AVERAGE ACROSS PERIODS
dat.per= dat.t3 %>% group_by(species, site, per, gdd.bin) %>% summarise(in6.cper=mean(in6.cper), in5.cper=mean(in5.cper), in4.cper=mean(in4.cper), in3.cper=mean(in3.cper), in2.cper=mean(in2.cper), in1.cper=mean(in1.cper), GDDs_binned=mean(GDDs_binned) )
  
#order by period
dat.per$per= factor(dat.per$per,levels=c("initial","cold","med","warm"), ordered=TRUE)

g1= ggplot(data=dat.per) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in4.cper, color="in4")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(site~per, scales = "free_x")+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in4.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in4.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,300)+ylab("proportional composition")

#----------------------
#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\DevelopInd\\")
file<-paste("Composition_boulderensis.pdf" ,sep="", collapse=NULL)
pdf(file,height = 10, width = 10)
g1
dev.off()

#--------------------
#DEVELOPMENTAL INDEX
#Plot DI by ordinal date

dat.b= subset(dat, dat$species %in% c("Melanoplus boulderensis") ) 

# AVERAGE ACROSS PERIODS
dat.bper= dat.b %>% group_by(site, elev, per, date.bin) %>% summarise(DI=mean(DI), ordinal=mean(ordinal) )
#order by period
dat.bper$per= factor(dat.bper$per,levels=c("initial","cold","med","warm"), ordered=TRUE)

#plot by date bin
di.plot.date= ggplot(data=dat.bper, aes(x=ordinal, y = DI, color=per))+facet_grid(site~.) +theme_bw()+geom_line(size=2)+xlim(140,220)

#by GDD bin
dat.bper= dat.b %>% group_by(site, elev, per, gdd.bin) %>% summarise(DI=mean(DI), cdd_sum=mean(cdd_sum) )
#order by period
dat.bper$per= factor(dat.bper$per,levels=c("initial","cold","med","warm"), ordered=TRUE)

#Plot DI by GDD
dat.bper$cdd= dat.bper$cdd_sum
di.plot.gdd= ggplot(data=dat.bper, aes(x=cdd, y = DI, color=per))+facet_grid(site~.)+theme_bw()+geom_line(size=2)+xlim(0,250)

#------------------
#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\DevelopInd\\")
file<-paste("DevelopmentalIndex_boulderensis.pdf" ,sep="", collapse=NULL)
pdf(file,height = 10, width = 10)
grid.arrange(di.plot.date,di.plot.gdd,ncol=2 )
dev.off() 
