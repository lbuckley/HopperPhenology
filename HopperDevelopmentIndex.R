 #load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

library(lme4)
library(nlme)
library(grDevices)

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

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )
hop= read.csv("HopperData.csv")

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

#--------------------
#BIN?

##Bin by GDD
gdds= seq( min(dat$GDDs, na.rm=TRUE), max(dat$GDDs, na.rm=TRUE), 30)
dat$gdd.bin= cut(dat$GDD, breaks = gdds, labels=FALSE)

dat.gddbin = dat %>% group_by(species,site,year,gdd.bin) %>% summarise_each(funs(mean))

##Bin by date
dates= seq( min(dat$ordinal, na.rm=TRUE)-1, max(dat$ordinal, na.rm=TRUE)+1, 20)
dat$date.bin= cut(dat$ordinal, breaks = dates, labels=FALSE, include.highest=TRUE)

dat.datebin = dat %>% group_by(species,site,year,date.bin) %>% summarise_each(funs(mean))

#dat$timeper= cut(dat$year, breaks = c(1957,1961,2011), labels=FALSE)
#------------------------------------------------------
#Composition plot

cols= gray.colors(4)
cols=rev(cols)

#Calculate percent composition
dat.t= dat.gddbin
dat.t= dat.t %>% mutate(in6.per= in6/total,in5.per= in5/total,in4.per= in4/total,in3.per= in3/total,in2.per= in2/total,in1.per= in1/total) 
#cumulative percentage
dat.t= dat.t %>% mutate(in6.cper= in6.per+in5.per+in4.per+in3.per+in2.per+in1.per, 
                        in5.cper= in5.per+in4.per+in3.per+in2.per+in1.per,
                        in4.cper= in4.per+in3.per+in2.per+in1.per,
                        in3.cper= in3.per+in2.per+in1.per,
                        in2.cper= in2.per+in1.per,
                        in1.cper= in1.per)

dat.t1= dat.t[,c(1:4,41:46) ] 
#To long format, needed?
dat.t2 <- gather(dat.t1, stage, cper, in6.cper:in1.cper, factor_key=TRUE)

#-------

#PLOT
dat.t3= subset(dat.t1, dat.t1$site=="B1")
dat.t3$GDDs_binned= gdds[dat.t3$gdd.bin]

g1= ggplot(data=dat.t3) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(species~year)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)

#--------------------------------
#Just boulderensis

dat.t3= subset(dat.t1, dat.t1$species=="Melanoplus boulderensis")
dat.t3$GDDs_binned= gdds[dat.t3$gdd.bin]

g1= ggplot(data=dat.t3) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(site~year)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+xlim(0,300)
