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

#----------------
#order by season gdd
#calculate number of seasonal gdds
clim.seas= clim %>% group_by(Year,Site) %>% summarise(dd.seas= max(cdd_sum) )
clim.seas$sy= paste(clim.seas$Site, clim.seas$Year, sep="_")

#plot seasonal dd by year
ggplot(data=clim.seas,aes(x=Year, y = dd.seas, color=Site)) + geom_point()+geom_line()

#------------
#group by seasonal GDD
#drop sites without data
clim.seas2= clim.seas
clim.seas2= clim.seas[-which(clim.seas$Site %in% c("A1")),]
clim.seas2= clim.seas[-which(clim.seas$Site %in% c("A1","CHA")),]

clim.seas2= clim.seas2 %>% group_by(Year) %>% summarise(dd.seas= mean(dd.seas) )
clim.seas2= clim.seas2[order(clim.seas2$dd.seas),]

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

#DEVELOPMENTAL INDEX
#Plot DI by ordinal date
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=per))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+geom_line()

#Plot DI by GDD
dat$cdd= dat$cdd_sum
di.plot= ggplot(data=dat, aes(x=cdd, y = DI, color=per))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+xlim(0,600)+ylim(1,6)+geom_line() #+geom_smooth(method="loess",se=F)
#note xlim restricted

#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenology\\figures\\")
pdf("DIplot.pdf",height = 10, width = 10)
di.plot
dev.off()

#-----------------
#SURFACE PLOT

#Add seasonal cdd
dat$sy=paste(dat$site, dat$year,sep="_")
match1= match(dat$sy, clim.seas$sy)
dat$dd.seas= NA
dat$dd.seas= unlist(clim.seas[match1,"dd.seas"])

#--------

#Interpolate
dat1=na.omit(dat)

#split by species and elevation
dat.mb= dat1[which(dat1$species=="Melanoplus boulderensis"),]

dat.e=  dat.mb[which(dat.mb$elev==2195),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.mb.2195 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.mb[which(dat.mb$elev==2591),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.mb.2591 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.mb[which(dat.mb$elev==3048),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.mb.3048 <- interp2xyz(s, data.frame=TRUE)

#--------------------------
dat.cp= dat1[which(dat1$species=="Camnula pellucida"),]

dat.e=  dat.cp[which(dat.cp$elev==2195),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.cp.2195 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.cp[which(dat.cp$elev==2591),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.cp.2591 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.cp[which(dat.cp$elev==3048),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.cp.3048 <- interp2xyz(s, data.frame=TRUE)

#--------------------------
dat.ms= dat1[which(dat1$species=="Melanoplus sanguinipes"),]

dat.e=  dat.ms[which(dat.ms$elev==1752),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.ms.1752 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.ms[which(dat.ms$elev==2195),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.ms.2195 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.ms[which(dat.ms$elev==2591),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.ms.2591 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.ms[which(dat.ms$elev==3048),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.ms.3048 <- interp2xyz(s, data.frame=TRUE)

#--------------------------
dat.md= dat1[which(dat1$species=="Melanoplus dawsoni"),]

dat.e=  dat.md[which(dat.md$elev==1752),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.md.1752 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.md[which(dat.md$elev==2195),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.md.2195 <- interp2xyz(s, data.frame=TRUE)

dat.e=  dat.md[which(dat.md$elev==2591),]
s=interp(x=dat.e$cdd,y=dat.e$dd.seas,z=dat.e$DI, duplicate="mean")
gdat.md.2591 <- interp2xyz(s, data.frame=TRUE)

#--------------------------
#plot
#mb
plot.di.mb.2195= ggplot(gdat.mb.2195) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "mb 2195m", size=5)

plot.di.mb.2591= ggplot(gdat.mb.2591) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "mb 2591", size=5)

plot.di.mb.3048= ggplot(gdat.mb.3048) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=200, label= "mb 3048", size=5)

#------------
#cp
plot.di.cp.2195= ggplot(gdat.cp.2195) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "cp 2195m", size=5)

plot.di.cp.2591= ggplot(gdat.cp.2591) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "cp 2591m", size=5)

plot.di.cp.3048= ggplot(gdat.cp.3048) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=200, label= "cp 3048m", size=5)

#------------
#ms
plot.di.ms.1752= ggplot(gdat.ms.1752) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=1000, label= "ms 1752m", size=5)

plot.di.ms.2195= ggplot(gdat.ms.2195) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "ms 2195m", size=5)

plot.di.ms.2591= ggplot(gdat.ms.2591) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "ms 2591m", size=5)

plot.di.ms.3048= ggplot(gdat.ms.3048) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=200, label= "ms 3048m", size=5)

#------------
#md
plot.di.md.1752= ggplot(gdat.md.1752) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=1000, label= "md 1752m", size=5)

plot.di.md.2195= ggplot(gdat.md.2195) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "md 2195m", size=5)

plot.di.md.2591= ggplot(gdat.md.2591) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="DI") +
  theme_bw(base_size=16)+xlab("gdd")+ylab("seasonal gdd")+theme(legend.position="none")+annotate("text", x=200,y=500, label= "md 2591m", size=5)

#-------------------
#PLOT
blank <- grid.rect(gp=gpar(col="white"))

di.plot <- grid.arrange(blank, plot.di.mb.2195, plot.di.mb.2591, plot.di.mb.3048, blank,plot.di.cp.2195, plot.di.cp.2591, plot.di.cp.3048,plot.di.ms.1752,plot.di.ms.2195, plot.di.ms.2591, plot.di.ms.3048,plot.di.md.1752,plot.di.md.2195, plot.di.md.2591, nrow=4)

#PLOT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenology\\figures\\")
pdf("DIplot.pdf",height = 12, width = 12)
plot(di.plot)
dev.off()

#==========================================================================
#Composition plot

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

#match
dat.t1$sy= paste(dat.t1$site, dat.t1$year, sep="_")
match1= match(dat.t1$sy,clim.seas$sy)
matched= which(!is.na(match1))

dat.t1$dd.seas=NA
dat.t1[matched,"dd.seas"]= clim.seas[match1[matched], "dd.seas"]
#round
dat.t1$dd.seas= round(dat.t1$dd.seas)

#------
#CODE BY GDD #see estimates below
#initial 1958,1959,1960
#cold 2009,2010,2014
#med 2008,2011,2013,2015
#warm 2006, 2007,2012

dat.t1$per=NA
dat.t1$per[which(dat.t1$year %in% c(1958,1959,1960) )]="initial"
dat.t1$per[which(dat.t1$year %in% c(2009,2010,2014) )]="cold"
dat.t1$per[which(dat.t1$year %in% c(2008,2011,2013,2015) )]="med"
dat.t1$per[which(dat.t1$year %in% c(2006,2007,2012) )]="warm"

#-------

#PLOT
dat.t3= subset(dat.t1, dat.t1$site=="B1")
dat.t3$GDDs_binned= gdds[dat.t3$gdd.bin]

g1= ggplot(data=dat.t3) + geom_line(aes(x=GDDs_binned, y = in5.cper, color="in5")) + geom_line(aes(x=GDDs_binned, y = in3.cper, color="in3")) +geom_line(aes(x=GDDs_binned, y = in2.cper, color="in2"))+facet_grid(species~year)+ geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + geom_ribbon(aes(x = GDDs_binned, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ geom_ribbon(aes(x = GDDs_binned, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)
