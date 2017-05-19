#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

sites= c("Redfox", "A1", "B1", "C1", "D1")  
elevs= c(1574, 2195, 2591, 3048, 3739)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

#--------------------------------------
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#==========================================================
#CODE FROM Buckley et al. 2014 PRSb

### SOME CODE CUT OUT IN THIS VERSION 
#BACK TO USING GDD BASED ON AIR TEMP
#account for snow melt

#see GDD binned plots
library(lme4)
library(nlme)
library(grDevices)

count= function(x){
  length( !is.na(x))
}

#load GDD data
#setwd("F:\\Work\\Grasshoppers\\JuvDev\\data\\GDD\\")
#Update date
setwd("F:\\Work\\Grasshoppers\\JuvDev\\data\\GDD\\GDDto2011\\")

A1gdd=read.csv("A1GDD_cleaned.csv")
B1gdd=read.csv("B1GDD.csv")
C1gdd=read.csv("C1GDD.csv")
D1gdd=read.csv("D1GDD.csv")
Chautgdd=read.csv("ChautauquaGDD.csv")

#load species data
setwd("F:\\Work\\Grasshoppers\\JuvDev\\data\\Processed\\Cleaned\\")

#Chaut
chaut.dat= read.csv("Chautdata.csv")
specs= c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus dawsoni", "Melanoplus sanguinipes")

#A1 
a1.dat= read.csv("A1data.csv")
specs= c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus dawsoni", "Melanoplus sanguinipes")
#also dodgei

#Fix NAs
inds=which(a1.dat$year=="1959" & is.na(a1.dat$GDD))
a1.dat[inds, "GDD"]= A1gdd[a1.dat[inds, "OrdinalDate"], "X1959"]
inds=which(a1.dat$year=="1960" & is.na(a1.dat$GDD))
a1.dat[inds, "GDD"]= A1gdd[a1.dat[inds, "OrdinalDate"], "X1960"]

#B1 
b1.dat= read.csv("B1data.csv")
ditch.dat= which(names(b1.dat)=="Genu_spec")
b1.dat= b1.dat[,-ditch.dat]
specs= c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus dawsoni", "Melanoplus dodgei")

#C1
c1.dat= read.csv("C1data.csv") 

#--------------------------------------
#UPDATE GDD
#CHAUT
gdd=Chautgdd
dat1= chaut.dat
na.ind= which(is.na(dat1$GDD))
dat.nagdd= dat1[na.ind,]
#match year
yrs= as.numeric(sub("X","",colnames(gdd)[3:ncol(gdd)]))
match.yr=match(dat.nagdd$year, yrs)+2

for(i in 1:nrow(dat.nagdd)){
  gdd1= gdd[ dat1[na.ind[i],"OrdinalDate"], match.yr[i] ]
  if(length(gdd1)>0) dat1[na.ind[i],"GDD"]= gdd1
}
chaut.dat=dat1

#A1
gdd=A1gdd
dat1= a1.dat
na.ind= which(is.na(dat1$GDD))
dat.nagdd= dat1[na.ind,]
#match year
yrs= as.numeric(sub("X","",colnames(gdd)[3:ncol(gdd)]))
match.yr=match(dat.nagdd$year, yrs)+2

for(i in 1:nrow(dat.nagdd)){
  gdd1= gdd[ dat1[na.ind[i],"OrdinalDate"], match.yr[i] ]
  if(length(gdd1)>0) dat1[na.ind[i],"GDD"]= gdd1
}
a1.dat=dat1

#B1
gdd=B1gdd
dat1= b1.dat
na.ind= which(is.na(dat1$GDD))
dat.nagdd= dat1[na.ind,]
#match year
yrs= as.numeric(sub("X","",colnames(gdd)[3:ncol(gdd)]))
match.yr=match(dat.nagdd$year, yrs)+2

for(i in 1:nrow(dat.nagdd)){
  gdd1= gdd[ dat1[na.ind[i],"OrdinalDate"], match.yr[i] ]
  if(length(gdd1)>0) dat1[na.ind[i],"GDD"]= gdd1
}
b1.dat=dat1

#C1
gdd=C1gdd
dat1= c1.dat
na.ind= which(is.na(dat1$GDD))
dat.nagdd= dat1[na.ind,]
#match year
yrs= as.numeric(sub("X","",colnames(gdd)[3:ncol(gdd)]))
match.yr=match(dat.nagdd$year, yrs)+2

for(i in 1:nrow(dat.nagdd)){
  gdd1= gdd[ dat1[na.ind[i],"OrdinalDate"], match.yr[i] ]
  if(length(gdd1)>0) dat1[na.ind[i],"GDD"]= gdd1
}
c1.dat= dat1

#-------------------------------------

#combine
dat= rbind(chaut.dat, a1.dat, b1.dat, c1.dat)

#get rid of entries with low sample size: at least three individuals, but keep if all first instar or adult
drop= which(dat$total<3 & dat$DI!=1 & dat$DI!=6)
dat=dat[-drop,]

#specs=unique(dat$Species)
years=sort(unique(dat$year))

#write.csv(dat,"dat_all.csv")

#add elevation
sites= c("Chaut", "A1", "B1", "C1")  
elevs= c(1752, 2195, 2591, 3048)

match1= match(dat$Site, sites)
dat$elev= elevs[match1]

#load length data
setwd("F:\\Work\\Grasshoppers\\JuvDev\\data\\")
lengths= read.csv("JuvenilesLengths_10Dec2011.csv")

dat.all=dat

#------------------------------------------------------------
#ANALYSIS

specs= c("Aeropedellus clavatus", "Camnula pellucida", "Melanoplus dawsoni", "Melanoplus sanguinipes", "Melanoplus dodgei")

#reduce to focal species
dat= dat[dat$Species %in% specs[2:5],]

#code period
dat$per=1
dat$per[dat$year>2000]=2

#code early and late species
dat$early_late=2
dat$early_late[dat$Species %in% specs[c(2,4)]]=1

#In relation to OD
mod1= lme(DI ~ poly(OrdinalDate,3)+elev*per* Species, random=~1|year, data=dat)
anova(mod1)

#GDD TEXT MODEL
dat=na.omit(dat)
#nlme model
mod1= lme(DI ~ poly(GDD,3)+elev*per+Species, random=~1|year, data=dat)
mod1= lme(DI ~ poly(GDD,3)+elev*per*Species, random=~1|year, data=dat)
anova(mod1)

dat1= subset(dat, dat$Species==specs[[2]])
mod1= lme(DI ~ poly(GDD,3)+elev*per, random=~1|year, data=dat1)

anova(mod1,type="marginal")

#------------------------------------------------------------
#PLOTS
#1. DI by GDD; DAYS BY STAGE
#2. COMP PLOT

#Plot development index across GDDs
# Past, Present; Binn GDDS

#####FIGURE 1
elevs_labs= paste(elevs,"m", sep="")

#BIN BY GDD AND ORDINAL DATE AND AVERAGE

#PLOT
#By GDD

setwd("F:\\Work\\HopperDevel\\figures\\Field\\")
file<-paste("DIbinned_ByGDD.pdf" ,sep="", collapse=NULL)
pdf(file,height = 8, width = 11)

par(mfrow=c(2,2), cex=1.2, mar=c(1.5, 2, 2, 0), mgp=c(2, 1, 0), oma=c(2,2,0,0), lwd=2)
#cols= c("gray90", "gray70", "gray50", "black") 
#colfunc <- colorRampPalette(c("green", "blue"))
#cols=colfunc(4)
cols=c("orange","limegreen","darkgreen","purple") 

dat= dat.all

#-------------------------
#clean data

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#Bin by GDD
gdds= seq( min(dat$GDD, na.rm=TRUE), max(dat$GDD, na.rm=TRUE), 30)
dat$gdd.bin= cut(dat$GDD, breaks = gdds, labels=FALSE)
dat$timeper= cut(dat$year, breaks = c(1957,1961,2011), labels=FALSE)
dat.gddagg.mean= aggregate(dat, by=list(dat$gdd.bin, dat$timeper, dat$Species, dat$Site), FUN='mean')
dat.gddagg.sd= aggregate(dat, by=list(dat$gdd.bin, dat$timeper, dat$Species, dat$Site), FUN='sd') 
dat.gddagg.count= aggregate(dat, by=list(dat$gdd.bin, dat$timeper, dat$Species, dat$Site), FUN='count')
names(dat.gddagg.mean)[1:4]= c("GDD.bin","TimePer", "Species", "Site")

dat.gddmat=data.frame(dat.gddagg.mean$Species, dat.gddagg.mean$Site, dat.gddagg.mean$gdd.bin,dat.gddagg.mean$timeper, dat.gddagg.mean$DI, dat.gddagg.sd$DI, dat.gddagg.count$DI)
colnames(dat.gddmat)=c("Species", "Site", "GDDbin","timeper", "mean", "sd", "N")
dat.gddmat = as.data.frame(dat.gddmat)
dat.gddmat$se= dat.gddmat$sd/sqrt(dat.gddmat$N)

#Remove incomplete data
#Remove M. sang B1
rem= which(dat.gddmat$Species==specs[4] & dat.gddmat$Site==sites[3])
dat.gddmat= dat.gddmat[-rem,]
#Remove M. daws A1
rem= which(dat.gddmat$Species==specs[3] & dat.gddmat$Site==sites[2])
dat.gddmat= dat.gddmat[-rem,]

#Bin by date
dates= seq( min(dat$OrdinalDate, na.rm=TRUE)-1, max(dat$OrdinalDate, na.rm=TRUE)+1, 20)
dat$date.bin= cut(dat$OrdinalDate, breaks = dates, labels=FALSE, include.highest=TRUE)
dat$timeper= cut(dat$year, breaks = c(1957,1961,2011), labels=FALSE)
dat.dagg.mean= aggregate(dat, by=list(dat$date.bin, dat$timeper, dat$Species, dat$Site), FUN='mean')
dat.dagg.sd= aggregate(dat, by=list(dat$date.bin, dat$timeper, dat$Species, dat$Site), FUN='sd') 
dat.dagg.count= aggregate(dat, by=list(dat$date.bin, dat$timeper, dat$Species, dat$Site), FUN='count')
names(dat.dagg.mean)[1:4]= c("date.bin","TimePer", "Species", "Site")  

dat.dmat=data.frame(dat.dagg.mean$Species, dat.dagg.mean$Site, dat.dagg.mean$date.bin,dat.dagg.mean$timeper, dat.dagg.mean$DI, dat.dagg.sd$DI, dat.dagg.count$DI)
colnames(dat.dmat)=c("Species", "Site","ODbin","timeper", "mean", "sd", "N")
dat.dmat = as.data.frame(dat.dmat)
dat.dmat$se= dat.dmat$sd/sqrt(dat.dmat$N)

for(sp in c(5,4,2,3)){ #dodg, sang, pell, daws
  
  #Plot by GDD
  dat.sp=subset(dat.gddmat, dat.gddmat$Species==specs[sp]) 
  
  for(si in 1: length(sites) ){    
    dat.si=subset(dat.sp, dat.sp$Site==sites[si])
    
    dat1=subset(dat.si, dat.si$timeper==1)
    dat2=subset(dat.si, dat.si$timeper==2)
    
    nas= which(is.na(dat1$Species))
    if(length(nas>0)) dat1= dat1[-nas,]
    nas= which(is.na(dat2$Species))
    if(length(nas>0)) dat2= dat2[-nas,]
    
    dat1= dat1[order(dat1$GDDbin),]
    dat2= dat2[order(dat2$GDDbin),]
    
    if(si==1) plot(gdds[dat1$GDDbin], dat1$mean, col=cols[si], ylab="DI", xlab="GDD", type="p",ylim=c(1,6), xlim=range(0,600),main=specs[sp], pch=16) 
    if(si>1) points(gdds[dat1$GDDbin], dat1$mean, col=cols[si], type="p", pch=16) 
    
    points(gdds[dat2$GDDbin], dat2$mean, col=cols[si],lty="dashed", type="p")
    
    #add smooth
    #LOESS
    # 	if(sum(is.na(gdds[dat1$GDDbin]))==0)try( points(loess.smooth(gdds[dat1$GDDbin], dat1$mean, family="gaussian"), type="l", col=cols[si]))
    #      if(sum(is.na(gdds[dat2$GDDbin]))==0)try( points(loess.smooth(gdds[dat2$GDDbin], dat2$mean, family="gaussian"), type="l", col=cols[si],lty="dashed"))
    # smooth.Pspline
    #spline, kernal
    
    ## Smoothing Spline
    #fit.sp = smooth.spline(dat1$mean ~ gdds[dat1$GDDbin], nknots=15)
    #lines(fit.sp, col="blue", lty=6)
    x1= gdds[dat1$GDDbin]
    y1= dat1$mean
    x2= gdds[dat2$GDDbin]
    y2= dat2$mean
    
    #set up for low sample sizes
    nx1= length(dat1$mean)
    nx2= length(dat2$mean)
    
    if(nx1>4 & sum(is.na(gdds[dat1$GDDbin]))==0 )( lines(smooth.spline(y1 ~ x1, nknots=min(nx1,8) ), col=cols[si])) #try
    if(nx2>4 & sum(is.na(gdds[dat2$GDDbin]))==0 )( lines(smooth.spline(y2 ~ x2, nknots=min(nx2,8) ), col=cols[si],lty="dashed")) #try
    
    if(nx1<=4) points(gdds[dat1$GDDbin], dat1$mean, col=cols[si], type="b", pch=16) 
    if(nx2<=4) points(gdds[dat2$GDDbin], dat2$mean, col=cols[si],lty="dashed", type="b")
    
    ## Natural Spline
    #library(splines)
    #fit.ns.3 = lm(y ~ ns(x, 3) )
    #lines(x, predict(fit.ns.3, data.frame(x=x)), col="green", lty=5)
    ## Kernel Curve
    #lines(ksmooth(x, y, "normal", bandwidth = 5), col = "orange", lty=7)
    
    #++++++++++++++++++++++++++
    
    #add se
    arrows(gdds[dat1$GDDbin], dat1$mean-dat1$se, gdds[dat1$GDDbin], dat1$mean+dat1$se, code=3, angle=90,length=0.1, col=cols[si])
    arrows(gdds[dat2$GDDbin], dat2$mean-dat2$se, gdds[dat2$GDDbin], dat2$mean+dat2$se, code=3, angle=90,length=0.1, col=cols[si])
    
    # add legend
    if(sp==5) legend(250,4, elevs_labs, col=cols, lty="solid", bty="n")
    
    # labels
    mtext("DI", side = 2, line = 0.5, outer = TRUE, cex=1.5)
    mtext("GDD", side = 1, line = 0.5, outer = TRUE, cex=1.5)
    
  }#end site loop
}#end species loop 

dev.off()

#------------------------------------------------------

############ FIGURE 2
#PLOT PERCENT COMPOSITION
cols= gray.colors(4)
cols=rev(cols)

#Calculate percent composition
#Sum across time periods within gdd.bin

#bin by date or GDD (date.bin)
dat.t= aggregate(dat[,c(5:10,12)], by=list(dat$timeper, dat$Species, dat$Site, dat$gdd.bin), FUN='sum', na.rm = TRUE)
colnames(dat.t)[1:4]=c("timeper", "Species", "Site", "date.bin")

#calculate GDD and date
dat.tm= aggregate(dat[,c(3,11)], by=list(dat$timeper, dat$Species, dat$Site, dat$gdd.bin), FUN='mean', na.rm = TRUE)
dat.t$OrdinalDate= dat.tm$OrdinalDate
dat.t$GDD= dat.tm$GDD

dat.percent= dat.t[,c("Stage1", "Stage2","Stage3","Stage4","Stage5","Stage6")]/dat.t$total
colnames(dat.percent)=c("PerStage1", "PerStage2","PerStage3","PerStage4","PerStage5","PerStage6")
dat.t= cbind(dat.t, dat.percent)

##Cummulative percentages
dat.t$cp.s1= dat.percent$PerStage1+dat.percent$PerStage2
#dat.t$cp.s2= dat.percent$PerStage1+dat.percent$PerStage2
dat.t$cp.s3= dat.percent$PerStage1+dat.percent$PerStage2+dat.percent$PerStage3+dat.percent$PerStage4
#dat.t$cp.s4= dat.percent$PerStage1+dat.percent$PerStage2+dat.percent$PerStage3+dat.percent$PerStage4
dat.t$cp.s5= dat.percent$PerStage1+dat.percent$PerStage2+dat.percent$PerStage3+dat.percent$PerStage4+dat.percent$PerStage5
dat.t$cp.s6= dat.percent$PerStage1+dat.percent$PerStage2+dat.percent$PerStage3+dat.percent$PerStage4+dat.percent$PerStage5+dat.percent$PerStage6

file<-paste("CompositionPlotFig.pdf" ,sep="", collapse=NULL)
pdf(file,height = 11, width = 11)
par(mfcol=c(4,2), cex=1.2, mar=c(2, 2, 2, 2), mgp=c(2, 1, 0), oma=c(2,2,0,0), lwd=2)

for(si in c(3,4) ){    #length(sites) #Compare B and C
  dat.si=subset(dat.t, dat.t$Site==sites[si])
  
  for(sp in c(5,2) ){ #loop species #C. pellucida and M. dodgei
    dat.sp=subset(dat.si, dat.si$Species==specs[sp])
    
    for(tp in 1:2){ #loop time periods
      dat2=subset(dat.sp, dat.sp$timeper==tp) 
      
      if(nrow(dat2)>0){ #check data exists
        
        #plot by ordinal date
        #dat2=dat2[order(dat2$OrdinalDate),]
        #x=dat2$OrdinalDate
        
        #plot by gdd
        dat2=dat2[order(dat2$GDD),]
        x=dat2$GDD
        
        times=c("1958-1960","2006-2011")
        spyr=paste(dat2$Species[1], times[dat2$timeper[1]], elevs_labs[si], sep=" ", collapse=NULL)
        plot(x, dat2$cp.s1, type="l", ylim=range(0,1), col="white", ylab="", xlab="", main=spyr, xlim=c(30,300))
        
        #shade 1st instaar
        xx <- c(x, rev(x) )
        yy <- c(dat2$cp.s1, rep(0,length(x)) )
        polygon(xx, yy, col=cols[1], border = "black")
        
        inds= match( c("cp.s1","cp.s3","cp.s5","cp.s6"), colnames(dat2) ) 
        #inds= match( c("cp.s1","cp.s2","cp.s3","cp.s4","cp.s5","cp.s6"), colnames(dat2) ) 
        
        for(i in 1:3){
          #shade quantiles
          xx <- c(x, rev(x))
          yy <- c(dat2[,inds[i]], rev(dat2[,inds[i+1]]) )
          polygon(xx, yy, col=cols[i+1], border = "black")
        }
        
      } #end check data exists
      
      #add legend
      ins.labs= c("1+2 Instars", "3+4 Instars", "5 Instar", "Adult" )
      if(si==4 & sp==5 & tp==1) legend(y=0.9,x=220, ins.labs, fill=cols, cex=0.8, ncol=1,bty="n")
      
    }}} #end timer per, site, species

mtext("Proportional Composition", side = 2, line = 0.5, outer = TRUE, cex=1.5)
mtext("GDD", side = 1, line = 0.5, outer = TRUE, cex=1.5)

dev.off()

#=============================================================
#Combine GDDS

#Combine data
cols= c("OrdinalDate","Site", "X1959", "X1960", "X2006", "X2007", "X2008") ##"X2009", "X2010", "X2011"
A1gdd$Site="A1"
B1gdd$Site="B1"
C1gdd$Site="C1"
D1gdd$Site="D1"

gdd.all=rbind(A1gdd[,cols], B1gdd[,cols], C1gdd[,cols], D1gdd[,cols])

#Col means
gdd.all$GDD.tp1= rowMeans( gdd.all[,c("X1959", "X1960")], na.rm = TRUE)
gdd.all$GDD.tp2= rowMeans( gdd.all[,c("X2006", "X2007", "X2008")], na.rm = TRUE)

#-------------------------
#FIGUE S1. PLOT GDD BY DATE 
#Ordinal dates May 15 to Sep 15, 135 to 258
par(mfrow=c(3,1), cex=1.2, mar=c(2, 2, 2, 1), mgp=c(2, 1, 0), oma=c(2,2,0,0), lwd=2)

#loop sites
for(si in 2:4){    
  gdd.si=subset(gdd.all, gdd.all$Site==sites[si])
  
  plot(gdd.si$OrdinalDate, gdd.si$GDD.tp2, lty="dashed", type="l", xlim=range(135,258),xlab="",ylab="", main=elevs_labs[si])
  points(gdd.si$OrdinalDate, gdd.si$GDD.tp1, type="l")
  
  if(si==2)legend(y=1000,x=130, c("1959-1960","2006-2011"), lty=c("solid","dashed"), cex=0.8, ncol=1,bty="n")
} #end loop sites

mtext("GDD", side = 2, line = 0.5, outer = TRUE, cex=1.5)
mtext("Ordinal Date", side = 1, line = 0.5, outer = TRUE, cex=1.5)

#---------------------------------------
#Figure 1 by ordinal date

#Bin by date
dates= seq( min(dat$OrdinalDate, na.rm=TRUE)-1, max(dat$OrdinalDate, na.rm=TRUE)+1, 20)
dat$date.bin= cut(dat$OrdinalDate, breaks = dates, labels=FALSE, include.highest=TRUE)
dat$timeper= cut(dat$year, breaks = c(1957,1961,2011), labels=FALSE)
dat.dagg.mean= aggregate(dat, by=list(dat$date.bin, dat$timeper, dat$Species, dat$Site), FUN='mean')
dat.dagg.sd= aggregate(dat, by=list(dat$date.bin, dat$timeper, dat$Species, dat$Site), FUN='sd') 
dat.dagg.count= aggregate(dat, by=list(dat$date.bin, dat$timeper, dat$Species, dat$Site), FUN='count')
names(dat.dagg.mean)[1:4]= c("ODbin","TimePer", "Species", "Site")  

dat.dmat=data.frame(dat.dagg.mean$Species, dat.dagg.mean$Site, dat.dagg.mean$date.bin,dat.dagg.mean$timeper, dat.dagg.mean$DI, dat.dagg.sd$DI, dat.dagg.count$DI)
colnames(dat.dmat)=c("Species", "Site","date.bin","timeper", "mean", "sd", "N")
dat.dmat = as.data.frame(dat.dmat)
dat.dmat$se= dat.dmat$sd/sqrt(dat.dmat$N)

setwd("F:\\Work\\HopperDevel\\figures\\Field\\")
file<-paste("DIbinned_ByOD.pdf" ,sep="", collapse=NULL)
pdf(file,height = 8, width = 11)

par(mfrow=c(2,2), cex=1.2, mar=c(1.5, 2, 2, 0), mgp=c(2, 1, 0), oma=c(2,2,0,0), lwd=2)

for(sp in c(5,4,2,3)){ #dodg, sang, pell, daws
  
  #Plot by OD
  dat.sp=subset(dat.dmat, dat.dmat$Species==specs[sp]) 
  
  for(si in 1: length(sites) ){    
    dat.si=subset(dat.sp, dat.sp$Site==sites[si])
    
    dat1=subset(dat.si, dat.si$timeper==1)
    dat2=subset(dat.si, dat.si$timeper==2)
    
    nas= which(is.na(dat1$Species))
    if(length(nas>0)) dat1= dat1[-nas,]
    nas= which(is.na(dat2$Species))
    if(length(nas>0)) dat2= dat2[-nas,]
    
    dat1= dat1[order(dat1$date.bin),]
    dat2= dat2[order(dat2$date.bin),]
    
    if(si==1) plot(dates[dat1$date.bin], dat1$mean, col=cols[si], ylab="DI", xlab="OD", type="p",ylim=c(1,6), xlim=range(100,260),main=specs[sp], pch=16) 
    if(si>1) points(dates[dat1$date.bin], dat1$mean, col=cols[si], type="p", pch=16) 
    
    points(dates[dat2$date.bin], dat2$mean, col=cols[si],lty="dashed", type="p")
    
    #add smooth
    
    ## Smoothing Spline
    x1= dates[dat1$date.bin]
    y1= dat1$mean
    x2= dates[dat2$date.bin]
    y2= dat2$mean
    
    #set up for low sample sizes
    nx1= length(dat1$mean)
    nx2= length(dat2$mean)
    
    if(nx1>4 & sum(is.na(dates[dat1$date.bin]))==0 )( lines(smooth.spline(y1 ~ x1, nknots=min(nx1,8) ), col=cols[si])) #try
    if(nx2>4 & sum(is.na(dates[dat2$date.bin]))==0 )( lines(smooth.spline(y2 ~ x2, nknots=min(nx2,8) ), col=cols[si],lty="dashed")) #try
    
    if(nx1<=4) points(dates[dat1$date.bin], dat1$mean, col=cols[si], type="b", pch=16) 
    if(nx2<=4) points(dates[dat2$date.bin], dat2$mean, col=cols[si],lty="dashed", type="b")
    
    #++++++++++++++++++++++++++
    
    #add se
    arrows(dates[dat1$date.bin], dat1$mean-dat1$se, dates[dat1$date.bin], dat1$mean+dat1$se, code=3, angle=90,length=0.1, col=cols[si])
    arrows(dates[dat2$date.bin], dat2$mean-dat2$se, dates[dat2$date.bin], dat2$mean+dat2$se, code=3, angle=90,length=0.1, col=cols[si])
    
    # add legend
    if(sp==5) legend(250,4, elevs_labs, col=cols, lty="solid", bty="n")
    
    # labels
    mtext("DI", side = 2, line = 0.5, outer = TRUE, cex=1.5)
    mtext("OD", side = 1, line = 0.5, outer = TRUE, cex=1.5)
    
  }#end site loop
}#end species loop 

dev.off()

#------------


