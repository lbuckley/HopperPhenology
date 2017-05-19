#load libraries
library(ggplot2)
library(plyr)
library(dplyr)

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
dat=na.omit(dat)
#nlme model
mod1= lme(DI ~ poly(GDDs,3)+elev*per+species, random=~1|year, data=dat)
mod1= lme(DI ~ poly(GDDs,3)+elev*per*species, random=~1|year, data=dat)
anova(mod1)

dat1= subset(dat, dat$species==specs[[2]])
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
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=3)+theme_bw()+geom_line()

#Plot DI by GDD
di.plot= ggplot(data=dat, aes(x=GDDs, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=1)+theme_bw()+geom_line()+xlim(0,600)
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
                        in2.cper= in2.per+in1.per)
                        
#PLOT

for(si in c(3,4) ){    #length(sites) #Compare B and C
  dat.si=subset(dat.t, dat.t$site==sites[si])
  
  for(sp in c(5,2) ){ #loop species #C. pellucida and M. dodgei
    dat.sp=subset(dat.si, dat.si$species==specs[sp])
    dat2=dat.sp
      
      if(nrow(dat2)>0){ #check data exists
        
        #plot by ordinal date
        #dat2=dat2[order(dat2$OrdinalDate),]
        #x=dat2$OrdinalDate
        
        #plot by gdd
        dat2=dat2[order(dat2$GDDs),]
        x=dat2$GDDs
        
        times=c("1958-1960","2006-2011")
        spyr=paste(dat2$species[1], elevs_labs[si], sep=" ", collapse=NULL)
        plot(x, dat2$in1.per, type="l", ylim=range(0,1), col="white", ylab="", xlab="", main=spyr, xlim=c(30,300))
        
        #shade 1st instaar
        xx <- c(x, rev(x) )
        yy <- c(dat2$in1.per, rep(0,length(x)) )
        polygon(xx, yy, col=cols[1], border = "black")
        
        inds= match( c("in1.per","in3.cper","in5.cper","in6.cper"), colnames(dat2) ) 
        
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
      
    }} #end site, species

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


