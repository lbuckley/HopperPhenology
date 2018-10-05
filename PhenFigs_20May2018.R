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
dat= dat[dat$species %in% specs,] #[3:6]
#order
dat$species= ordered(dat$species, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Chloealtis abdominalis", "Camnula pellucida","Melanoplus sanguinipes","Melanoplus dawsoni") )

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
dat$species= factor(dat$species, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Chloealtis abdominalis", "Camnula pellucida", "Melanoplus sanguinipes", "Melanoplus dawsoni") )

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
#==============================================
#FIGURE 1
#ELEVATION PLOTS

#plot
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenology/figures/")

#check years with GDD data
clim2= aggregate(clim1, list(clim1$Year), FUN=mean)

#subset to years with data for all
clim2= subset(clim1, clim$Year %in% c(1958:1969, 2006:2011,2014) )

plot_gdd_elev=ggplot(data=clim2, aes(x=elevation, y = Cdd_seas, group=Year, color=Cdd_siteave, linetype=period))+
  geom_line()+ #geom_point()+
  scale_colour_gradientn(name="mean season gdd", colours =matlab.like(10))+
  theme_bw()+ylab("season growing degree days (C)")+xlab("elevation (m)")

#-------
#Plot ave phenology across years

#elevation to numeric
dat$elevation= as.numeric(as.character(dat$elev))

#aggregte to spsiteyear
dat.ave= aggregate(dat, list(dat$spsiteyear, dat$species, dat$site, dat$year, dat$siteyear),FUN=mean, na.rm=TRUE)
names(dat.ave)[1:5]= c('spsiteyear', 'species', 'site', 'year','siteyear')

#check siteyear
dat.ave2= aggregate(dat.ave, list(dat.ave$siteyear), FUN=length)
#restrict to years with data across sites, 2007:2011
dat.ave= subset(dat.ave, dat.ave$year %in% 2007:2011)

#average across years
dat.ave= aggregate(dat.ave, list(dat.ave$elevation, dat.ave$species),FUN=mean, na.rm=TRUE )
dat.ave$species= dat.ave$Group.2

plot_gdds_elev=ggplot(data=dat.ave, aes(x=elevation, y = gdd_adult, group=species, color=species))+
  geom_line()+ geom_point()+
  theme_bw()+ylab("cummulative growing degree days (C)")+xlab("elevation (m)")

plot_doy_elev=ggplot(data=dat.ave, aes(x=elevation, y = doy_adult, group=species, color=species))+
  geom_line()+ geom_point()+
  theme_bw()+ylab("day of year")+xlab("elevation (m)")+theme(legend.position="none")

#check data
unique(dat$siteyear)

#----------

#plot together
fig0= plot_grid(plot_gdd_elev, plot_doy_elev, plot_gdds_elev, labels = c('A', 'B','C'), rel_widths=c(1.2,0.75,1.3), nrow=1)

pdf("Fig1_GDD_phen_byElev.pdf", height = 5, width = 12)
fig0
dev.off()

#====================================
## FIGURE 2.
 
#DEVELOPMENTAL INDEX
#Plot DI by ordinal date

#update elevation labels
dat$elev.lab= paste(dat$elev,"m",sep="")
dat$elev.lab= factor(dat$elev.lab, levels=c("3048m","2591m","2195m","1752m") )

di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=Cdd_siteave, group=siteyear, linetype=period))+facet_grid(elev.lab~species) +
  theme_bw()+
  geom_point()+geom_line(aes(alpha=0.5))+ #+geom_smooth(se=FALSE, aes(alpha=0.5), span=2)+
  scale_colour_gradientn(colours =matlab.like(10))+ylab("development index")+xlab("day of year")+labs(color="mean season gdds")+
  theme(legend.position = "bottom") + guides(alpha=FALSE)

#Plot DI by GDD
di.plot.gdd= ggplot(data=dat, aes(x=cdd_sum, y = DI, color=Cdd_siteave, group=siteyear, linetype=period))+facet_grid(elev.lab~species) +
  theme_bw()+
  geom_point()+geom_line(aes(alpha=0.5))+ #+geom_smooth(se=FALSE, aes(alpha=0.5),span=2)+
  scale_colour_gradientn(colours =matlab.like(10))+ylab("development index")+xlab("cummulative growing degree days")+labs(color="mean season gdds")+
  xlim(0,1100)+
  theme(legend.position = "bottom") + guides(alpha=FALSE)

#----
pdf("Fig2_plot_DI.pdf", height = 10, width = 12)
di.plot
dev.off()

pdf("FigSX_plot_DIgdd.pdf", height = 10, width = 12)
di.plot.gdd
dev.off()

#--------
#Calculate r^2 and slope of polynomial
#across years
#stop when reach adulthood

#combination of sp and sites
dat$spsite= paste(dat$species, dat$site, sep="_")
spsites= unique(dat$spsite)

rs= as.data.frame(matrix(NA, length(spsites), 5))
rs[,1]=spsites
rs[,4]=matrix(unlist(strsplit(spsites, split="_")),ncol=2, byrow=T)[,1]
rs[,5]=matrix(unlist(strsplit(spsites, split="_")),ncol=2, byrow=T)[,2]

for(ind in 1:length(spsites) ){
 dat.ss= subset(dat, dat$spsite==spsites[ind])
 dat.ss$yeardi= paste(dat.ss$year, as.character(round(dat.ss$DI,2)), sep="" )
#sort by ordinal date
 dat.ss=dat.ss[order(dat.ss$yeardi),]
 #drop data once reach DI=6
 dat.ss=dat.ss[ (duplicated(dat.ss$yeardi)==TRUE & dat.ss$DI>=6)==FALSE,]
 
 #calculate slope and r2
modr= summary(lm(DI~poly(cdd_sum,2,raw=TRUE), data=dat.ss ))
rs[ind,2]= modr$coefficients[2,1]
rs[ind,3]= modr$r.squared
}

#--------------------
#plot with variable x axis

di.plot.gdd.free= ggplot(data=dat, aes(x=cdd_sum, y = DI, color=Cdd_siteave, group=siteyear, linetype=period))+
  facet_grid(species~elev.lab, scales="free_x") +
  theme_bw()+
  geom_point()+geom_line(aes(alpha=0.5))+ #+geom_smooth(se=FALSE, aes(alpha=0.5),span=2)+
  scale_colour_gradientn(colours =matlab.like(10))+ylab("development index")+xlab("cummulative growing degree days")+labs(color="mean season gdds")+
  theme(legend.position = "bottom") + guides(alpha=FALSE)+ geom_text(aes(x=200, y=5, label=V2, group=NULL), color="black", data=rs) 
##FIX

#--------
pdf("FigSX_plot_DIgdd_freex.pdf", height = 10, width = 12)
di.plot.gdd.free
dev.off()

#====================================
## FIGURE 3.
# PLOT ADULT PHENOLOGY
# estimated by DI

#aggregate to spsiteyear
dat.ssy= dat[,c("elev","species","cdd_seas","doy_adult","gdd_adult","spsiteyear","elev.lab","elevation","period")  ]
dups= duplicated(dat.ssy$spsiteyear)
dat.ssy=dat.ssy[which(dups==FALSE),]

#PLOT
#calculate significant regressions
#apply through combinations of species and elevations
dat.ssy$elevspec= paste(dat.ssy$elev, dat.ssy$species, sep="")
elevspec= matrix(unique(dat.ssy$elevspec))

#---
#EXTRACT COEFFS
p.doy= apply(elevspec,1, FUN=function(x) summary(lm(dat.ssy$doy_adult[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)))] ~ dat.ssy$cdd_seas[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)) )]) )$coefficients[2,])
p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(dat.ssy$gdd_adult[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)))] ~ dat.ssy$cdd_seas[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)) )]) )$coefficients[2,])
#combine
p.mat=as.data.frame(cbind(elevspec,t(p.doy),t(p.gdd) ))
#make numeric
p.mat$Estimate= as.numeric(as.character(p.mat$Estimate))
p.mat$`Std. Error`= as.numeric(as.character(p.mat$`Std. Error`))
#CIs
p.mat$upCI= p.mat$Estimate + 1.96*p.mat$`Std. Error`
p.mat$lowCI= p.mat$Estimate - 1.96*p.mat$`Std. Error`
#---

#extract p-values
p.doy= apply(elevspec,1, FUN=function(x) summary(lm(dat.ssy$doy_adult[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)))] ~ dat.ssy$cdd_seas[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(dat.ssy$gdd_adult[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)))] ~ dat.ssy$cdd_seas[which(dat.ssy$elev==substr(x,1,4)&dat.ssy$species==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
#combine
p.mat=as.data.frame(cbind(elevspec,p.doy,p.gdd))
#add columns for significance
p.mat$sig.doy<-"nonsignificant"
p.mat$sig.doy[which(p.doy<0.05)]="significant"
p.mat$sig.gdd="nonsignificant"
p.mat$sig.gdd[which(p.gdd<0.05)]="significant"

#add back to matrix
match1= match(dat.ssy$elevspec, elevspec)
dat.ssy$sig.doy= factor(p.mat[match1,"sig.doy"], levels=c("significant","nonsignificant"))
dat.ssy$sig.gdd= factor(p.mat[match1,"sig.gdd"], levels=c("significant","nonsignificant"))
#---

#DOY
plot.phen.doye=ggplot(data=dat.ssy, aes(x=cdd_seas, y = doy_adult, color=species))+
  geom_point(aes(shape=period, fill=species, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, fill=NULL, stroke=1), size=3)+
  geom_smooth(method="lm",se=F, aes(linetype=sig.doy))+
  facet_wrap(~elev.lab, ncol=1, scales="free") +
  theme_bw()+ylab("day of year")+xlab("season growing degree days (C)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))+theme(legend.position="none")

#GDD
plot.phen.gdde=ggplot(data=dat.ssy, aes(x=cdd_seas, y = gdd_adult, color=species))+
  geom_point(aes(shape=period, fill=species, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, fill=NULL, stroke=1), size=3)+
  geom_smooth(method="lm",se=F, aes(linetype=sig.gdd))+
  facet_wrap(~elev.lab, ncol=1, scales="free") +
  theme_bw()+ylab("cummulative growing degree days")+xlab("season growing degree days (C)")+
  labs(linetype="significance")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))

pdf("Fig3_phen_est.pdf", height = 12, width = 10)
plot_grid(plot.phen.doye, plot.phen.gdde, nrow=1, rel_widths=c(1,1.5) )
dev.off()

#=================================
#COMPOSITION PLOT

#add elevation
dat$elevation= elevs[match(dat$site, sites)]

##Bin by GDD
dat$GDDs<- dat$cdd_sum
gdds= seq( min(dat$GDDs, na.rm=TRUE), max(dat$GDDs, na.rm=TRUE), length.out=20) #CHANGE NUMBER GDD BINS
dat$gdd.bin= cut(dat$GDDs, breaks = gdds, labels=FALSE)

dat.gddbin = dat %>% group_by(species,site,year,gdd.bin) %>% summarise_each(funs(mean))

##Bin by date
dates= seq( min(dat$ordinal, na.rm=TRUE)-1, max(dat$ordinal, na.rm=TRUE)+1, 20)
dat$date.bin= cut(dat$ordinal, breaks = dates, labels=FALSE, include.highest=TRUE)

dat.datebin = dat %>% group_by(species,site,year,date.bin) %>% summarise_each(funs(mean))
#----------------

#Calculate percent composition
dat.t= dat  #dat.datebin

dat.t= dat.t %>% mutate(in6.per= in6/total,in5.per= in5/total,in4.per= in4/total,in3.per= in3/total,in2.per= in2/total,in1.per= in1/total) 
#cumulative percentage
dat.t= dat.t %>% mutate(in6.cper= in6.per+in5.per+in4.per+in3.per+in2.per+in1.per, 
                        in5.cper= in5.per+in4.per+in3.per+in2.per+in1.per,
                        in4.cper= in4.per+in3.per+in2.per+in1.per,
                        in3.cper= in3.per+in2.per+in1.per,
                        in2.cper= in2.per+in1.per,
                        in1.cper= in1.per)

#-------
dat.t3<- dat.t

#PLOT
#dat.t3= subset(dat.t1, dat.t1$site %in% c("A1","C1"))
dat.t3= subset(dat.t3, dat.t3$year %in% c(2010, 2007))

#make warm and cold labels
dat.t3$tempyear<- "2010 cool"
dat.t3$tempyear[which(dat.t3$year==2007)]<-"2007 warm"
dat.t3$tempyear= factor(dat.t3$tempyear, levels=c("2010 cool","2007 warm"))

# #dat.t3$GDDs_binned= gdds[dat.t3$gdd.bin]
# dat.t3$GDDs_binned= dates[dat.t3$date.bin]
# 
# dat.t3$per= factor(dat.t3$per, levels=c("initial","cold","med","warm") )

dat.t3$species= factor(dat.t3$species, levels=c("Melanoplus boulderensis","Camnula pellucida","Melanoplus sanguinipes","Melanoplus dawsoni"))

#elevation factor
dat.t3$elevation= factor(dat.t3$elevation, levels=c(3048,2591,2195,1752) )
dat.t3$Cdd_siteave= factor(round(dat.t3$Cdd_siteave,2))

#restrict columns to remove NAs
dat.t3= dat.t3[,c("ordinal","species","elevation","Cdd_siteave","year","tempyear","in6.cper","in5.cper","in4.cper","in3.cper","in2.cper","in1.cper","cdd_sum")]

#SUBSET SPECIES
#dat.t3= subset(dat.t3, dat.t3$species %in% c("Melanoplus boulderensis")) #,"Camnula pellucida", "Melanoplus boulderensis", "Melanoplus dawsoni"

comp.md= ggplot(data=dat.t3[which(dat.t3$species=="Melanoplus dawsoni"),]) + geom_line(aes(x=ordinal, y = in5.cper, color="4th and 5th")) + 
  geom_line(aes(x=ordinal, y = in3.cper, color="3rd")) +
  geom_line(aes(x=ordinal, y = in2.cper, color="1st and 2nd"))+
  facet_grid(elevation~tempyear)+xlim(135,220)+
  geom_ribbon(aes(x = ordinal, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + 
  geom_ribbon(aes(x = ordinal, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + 
  geom_ribbon(aes(x = ordinal, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ 
  geom_ribbon(aes(x = ordinal, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+
  ylab("proportion composition")+xlab("day of year")+
  theme(legend.position="bottom", text = element_text(size=14))+labs(color="Instar" )+  
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_vline(xintercept = 152, color="gray")+geom_vline(xintercept = 213, color="gray")
#+scale_color_manual(values=c("red","green","purple"),labels=c("1st and 2nd","3rd","4th and 5th") )

comp.mb= ggplot(data=dat.t3[which(dat.t3$species=="Melanoplus boulderensis"),]) + geom_line(aes(x=ordinal, y = in5.cper, color="4th and 5th")) + 
  geom_line(aes(x=ordinal, y = in3.cper, color="3rd")) +
  geom_line(aes(x=ordinal, y = in2.cper, color="1st and 2nd"))+
  facet_grid(elevation~tempyear)+xlim(135,220)+
  geom_ribbon(aes(x = ordinal, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + 
  geom_ribbon(aes(x = ordinal, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + 
  geom_ribbon(aes(x = ordinal, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ 
  geom_ribbon(aes(x = ordinal, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+
  ylab("proportion composition")+xlab("day of year")+
  theme(legend.position="bottom", text = element_text(size=14))+labs(color="Instar" )+  
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_vline(xintercept = 152, color="gray")+geom_vline(xintercept = 213, color="gray")

## FIGURE SX. Composition
#PLOT
pdf("Fig4_composition.pdf",height = 8, width = 12)
plot_grid(comp.mb, comp.md, nrow=1, rel_widths=c(1,1),labels = c('A) M. boulderensis','B) M. dawsoni'))
dev.off()

#-------------------------
#COMPOSITION GDD VERSION

comp.md.gdd= ggplot(data=dat.t3[which(dat.t3$species=="Melanoplus dawsoni"),]) + geom_line(aes(x=cdd_sum, y = in5.cper, color="4th and 5th")) + 
  geom_line(aes(x=cdd_sum, y = in3.cper, color="3rd")) +
  geom_line(aes(x=cdd_sum, y = in2.cper, color="1st and 2nd"))+
  facet_grid(elevation~tempyear)+xlim(0,900)+
  geom_ribbon(aes(x = cdd_sum, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + 
  geom_ribbon(aes(x = cdd_sum, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + 
  geom_ribbon(aes(x = cdd_sum, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ 
  geom_ribbon(aes(x = cdd_sum, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+
  ylab("proportion composition")+xlab("cummulative growing degree days (C)")+
  theme(legend.position="bottom", text = element_text(size=14))+labs(color="Instar" )+  
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_vline(xintercept = 0, color="gray")+geom_vline(xintercept = 800, color="gray")
#+scale_color_manual(values=c("red","green","purple"),labels=c("1st and 2nd","3rd","4th and 5th") )

comp.mb.gdd= ggplot(data=dat.t3[which(dat.t3$species=="Melanoplus boulderensis"),]) + geom_line(aes(x=cdd_sum, y = in5.cper, color="4th and 5th")) + 
  geom_line(aes(x=cdd_sum, y = in3.cper, color="3rd")) +
  geom_line(aes(x=cdd_sum, y = in2.cper, color="1st and 2nd"))+
  facet_grid(elevation~tempyear)+xlim(0,900)+
  geom_ribbon(aes(x = cdd_sum, ymin = in2.cper, ymax =1),fill = "orange", alpha = 0.4) + 
  geom_ribbon(aes(x = cdd_sum, ymin = in3.cper, ymax = in5.cper),fill = "blue", alpha = 0.4) + 
  geom_ribbon(aes(x = cdd_sum, ymin = in2.cper, ymax = in3.cper),fill = "green", alpha = 0.4)+ 
  geom_ribbon(aes(x = cdd_sum, ymin = 0, ymax = in2.cper),fill = "red", alpha = 0.4)+
  ylab("proportion composition")+xlab("cummulative growing degree days (C)")+
  theme(legend.position="bottom", text = element_text(size=14))+labs(color="Instar" )+  
  guides(colour = guide_legend(override.aes = list(size=3)))+
  geom_vline(xintercept = 0, color="gray")+geom_vline(xintercept = 800, color="gray")

## FIGURE SX. Composition
#PLOT
pdf("Fig4_composition_gdd.pdf",height = 8, width = 12)
plot_grid(comp.mb.gdd, comp.md.gdd, nrow=1, rel_widths=c(1,1),labels = c('A) M. boulderensis','B) M. dawsoni'))
dev.off()

#----------------

#Duration

hop4$siteyearsp= paste(hop4$siteyear, hop4$species, sep="")
hop4$siteyearspq= paste(hop4$siteyear, hop4$species,hop4$quantile, sep="")
siyrsp= unique(hop4$siteyearsp)

siyrsp.up= paste(siyrsp, "85", sep="")
siyrsp.low= paste(siyrsp, "15", sep="")

#match upper and lower quantiles
match.up=match(siyrsp.up, hop4$siteyearspq) 
match.low=match(siyrsp.low, hop4$siteyearspq) 

#diff
hop.diff= data.frame(siyrsp, dur=hop4[match.up, "ordinal"]-hop4[match.low, "ordinal"])
match1= match(hop4$siteyearsp, siyrsp)
hop4$diff= hop.diff[match1,"ordinal"]

#plot duration
plot.dur=ggplot(data=hop4, aes(x=cdd_seas, y = diff, color=species))+
  geom_point(aes(shape=period, fill=species, alpha=period, stroke=1), size=3)+
  geom_point(aes(shape=period, fill=NULL, stroke=1), size=3)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~elevation, ncol=1, scales="free") +
  theme_bw()+ylab("duration of adulthood")+xlab("season growing degree days (C)")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_alpha_manual(values = c(0.2,0.9))

pdf("Fig_duration.pdf",height = 8, width = 5)
plot.dur
dev.off()

#--------------

