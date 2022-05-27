library(reshape)
library(reshape2)
library(tidyr)
library(cowplot)
require(nlme)
require(lme4)
require(plyr)
library(MuMIn)
library(RColorBrewer)
library(lemon)
library(ggplot2)

## DATA WITH ALL NYMPHS
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
dat.all= read.csv("HopperData_Sept2019_forPhenOverlap.csv")

#-----------------

#LOAD SPECIES DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
spec.diapause= read.csv("SpeciesDiapause.csv")
#dat.all1= read.csv("HopperClimateData.csv")
#data currently from phenology file PhenFigs_20May2018.R
timing.mat= read.csv("timingmat_2019.csv" )
colnames(timing.mat)[1]<-"species"
rownames(timing.mat)<- timing.mat$species
#alphabetize by group
timing.mat= timing.mat[order(timing.mat$timing, timing.mat$species),]

#-------------
#climate data
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/climate/")
#read annual precipitation and snow data (not site specific)
clim.ann= read.csv("AlexanderYearlyClimate.csv")

#read data for all sites
clim= read.csv("AlexanderClimateAll_filled_May2022.csv")
#yearly data
clim.ave= aggregate(clim, by=list(clim$Year, clim$Site), FUN="mean", na.rm=TRUE)

#------------------------------------
#compare egg diapausers vs non
dat.all$diapause="egg"
dat.all$species= as.character(dat.all$species)
dat.all$diapause[which(dat.all$species %in% c("Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes"))]="nymph"

#fix names
dat.all$species[which(dat.all$species=="Melanoplus bivatattus")]="Melanoplus bivittatus"

#focal species
specs= c("Aeropedellus clavatus","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes","Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes","Melanoplus fasciatus") #"Chloealtis abdominalis",
specs= unique(dat.all$species)

#subset to focal species
##dat= subset(dat.all, dat.all$species %in% specs)
#keep all species
dat=dat.all

years= unique(na.omit(dat$year))
sites= unique(dat$site)

#elevs and sites
sites=c("CHA", "A1", "B1", "C1")
elevs= c(1752, 2195, 2591, 3048)

#Calculate development index
dat$DIp=0
inds=which(dat$total>0)  
dat$DIp[inds]= (dat$in1[inds] +dat$in2[inds]*2 +dat$in3[inds]*3 +dat$in4[inds]*4 +dat$in5[inds]*5 +dat$in6[inds]*6)/6

dat$spsiteyear= paste(dat$siteyear, dat$species, sep="")

#calculate proportion
dat.sum= ddply(dat, c("site", "year","species","spsiteyear"), summarise,
               DIp = sum(DIp, na.rm=TRUE), TotMean=mean(total, na.rm=TRUE), TotSum=sum(total, na.rm=TRUE), AdultSum=sum(in6, na.rm=TRUE), Ndates=length(total[DIp>0]) )
dat.sum= dat.sum[order(dat.sum$site,dat.sum$species, dat.sum$year),]

dat$DIptotal=NA
dat$Ndates= NA
dat$TotSum= NA
dat$AdultSum= NA
match1= match(dat$spsiteyear, dat.sum$spsiteyear)
matched= which(!is.na(match1))
dat$DIptotal[matched]<- dat.sum$DIp[match1[matched]]  
dat$Ndates[matched]<- dat.sum$Ndates[match1[matched]]
dat$TotSum[matched]<- dat.sum$TotSum[match1[matched]]
dat$AdultSum[matched]<- dat.sum$AdultSum[match1[matched]]
#normalize
dat$DIpNorm= dat$DIp/dat$DIptotal
dat$AdultNorm= dat$in6/dat$AdultSum

#drop cases with TotSum<30 or Ndates<5
dat= dat[which(dat$Ndates>5),]
dat= dat[which(dat$TotSum>30),]

#========================================
#Fit distributions
#fit beta: mode, breadth, skewness

#https://stats.stackexchange.com/questions/222233/residual-standard-error-difference-between-optim-and-glm
#https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet

#fit gaussian
fitG =
  function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    fit= optim(c(mu,sig,scale),f, method="L-BFGS-B", lower=c(100,20,50), upper=c(300,200,5000), hessian=TRUE)
    estimates <- fit$par     # Parameters estimates
    se <- sqrt(diag(solve(fit$hessian))) # Standard errors of the estimates
    n=length(x)
    sd=se*n^2
    
    return(c(estimates,sd,fit$convergence))
    
      }
#fit= fitG(temp,y,mu=200, sig=1, scale=1)

#data storage
fits= array(data=NA, dim=c(length(specs), length(years), length(sites),7), list(specs,years,sites,c("mu","sig","scale","mu.sd","sig.sd","scale.sd","abund")) )

colors= rainbow(9)

setwd("/Volumes/GoogleDrive/Shared drives/RoL_FitnessConstraints/projects/Extremes/figures/")

pdf("Fits.pdf", height = 10, width = 10)
par(mfrow=c(5,5))

for(site in 1:length(sites)){
  dat.sub= dat[which(dat$site==sites[site] ),]
  
  spec.inds= which(specs %in% dat.sub$species )
  year.inds= which(years %in% dat.sub$year )
  
  for(spec.ind in spec.inds){

  for(year.ind in year.inds){
  
      dat.sub2= dat.sub[which(dat.sub$species==specs[spec.ind] & dat.sub$year==years[year.ind]),]
  
      plot= FALSE
      if(nrow(dat.sub2)>5){
      p= fitG(dat.sub2$ordinal, dat.sub2$DIp, mu=150, sig=50, scale=100 )
      if(p[7]==0) fits[spec.ind, year.ind, site, 1:6]= p[1:6]
      
      if(year.ind==1 | plot==FALSE) plot(dat.sub2$ordinal, dat.sub2$DIp, type="p", xlim=range(100,300), ylim=range(0,100), col= colors[year.ind], main= paste(sites[site],specs[spec.ind],sep=" ") )
      plot=TRUE
      if(year.ind>1) points(dat.sub2$ordinal, dat.sub2$DIp, type="p", col= colors[year.ind])
      #plot fit
      lines(dat.sub2$ordinal,p[3]*dnorm(dat.sub2$ordinal,p[1],p[2]), col= colors[year.ind])
      
      #add abundance
      fits[spec.ind, year.ind, site, 7]= dat.sub2$AdultSum[1]
      }
    }
  }
}
dev.off()

#----------------
#expand array
fit.df=as.data.frame.table(fits)
colnames(fit.df)=c("species","year","site","param","value")
#make numeric and add elevation
fit.df$year= as.numeric(as.character(fit.df$year))
fit.df$elev<- elevs[fit.df$site]
fit.df= fit.df[!is.na(fit.df$value),]
fit.df= fit.df[which(fit.df$value>0),]

#add period
fit.df$period<- "resurvey"
fit.df$period[which(fit.df$year<1965)]<- "initial"

#add diapause
fit.df$diapause="egg"
fit.df$diapause[which(fit.df$species %in% c("Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes"))]="nymph"

#add gdds
fit.df$siteyear= paste(fit.df$site, fit.df$year, sep="")
dat$siteyear= paste(dat$site, dat$year, sep="")

match1= match(fit.df$siteyear, dat$siteyear)
fit.df$Cdd_siteave= dat$Cdd_siteave[match1]
fit.df$cdd_seas= dat$cdd_seas[match1]
fit.df$Tmean= dat$Tmean[match1]

#add seasonal timing (calculated below)
fit.df$timing= timing.mat$timing[match(fit.df$species, timing.mat$species)]
#drop species without timing
fit.df= fit.df[which(!is.na(fit.df$timing)),]

#order species by timing
fit.df$species= factor(fit.df$species, levels=timing.mat$species )

#add elevation label
fit.df$elev.lab= paste(fit.df$elev, "m", sep="")

#add annual snow and precip data
match1= match(fit.df$year, clim.ann$Year)
fit.bind= clim.ann[match1, c("pre","snow","SpringPre","SummerPre","SpringSnow")]
fit.df= cbind(fit.df, fit.bind)

#add extremes
exts= rbind( cbind(rep("CHA", nrow(clim.ann)), clim.ann$Year, clim.ann$CHA_tx90, clim.ann$CHA_wsdi),
  cbind(rep("A1", nrow(clim.ann)), clim.ann$Year, clim.ann$A1_tx90, clim.ann$A1_wsdi),
  cbind(rep("B1", nrow(clim.ann)), clim.ann$Year, clim.ann$B1_tx90, clim.ann$B1_wsdi),
  cbind(rep("C1", nrow(clim.ann)), clim.ann$Year, clim.ann$C1_tx90, clim.ann$C1_wsdi) )
exts= as.data.frame(exts)
colnames(exts)=c("Site","Year","Tx90","WSDI")
exts$SiteYr= paste(exts$Site, exts$Year, sep="")

match1= match(fit.df$siteyear, exts$SiteYr)
fit.df$Tx90= as.numeric(exts$Tx90[match1])
fit.df$WSDI= as.numeric(exts$WSDI[match1])

#cast
rownames(fit.df) <- c()
fit.df1 <- dcast(fit.df, species + year + site + Cdd_siteave +Tmean ~ param, value.var="value")
fit.df1$spsiteyear= paste(fit.df1$species,fit.df1$year, fit.df1$site, sep="_")

#--------------------
# Compare environmental metrics to fits

fit.df$param.lab="peak of abundance distribution"
fit.df$param.lab[which(fit.df$param=="sig")]= "breadth of abundance distribution"
fit.df$param.lab[which(fit.df$param=="scale")]= "scale of abundance distribution"
fit.df$param.lab= factor(fit.df$param.lab, levels=c("peak of abundance distribution","breadth of abundance distribution","scale of abundance distribution") )

#set species colors
levels(fit.df$species)
cols= c(brewer.pal(4,"Blues")[2:4],  brewer.pal(5,"Greens")[4:5], brewer.pal(9,"YlOrRd") )
#image(1:14,1,as.matrix(1:14),col=cols)

fit.df$elev= factor(fit.df$elev)

#predictors
#year, pre, snow, Tmean, cdd_seas, SpringPre, SummerPre, SpringSnow, Tx90, WSDI

ggplot(fit.df[which(fit.df$param %in%c("mu","sig","scale","abund") & fit.df$species=="Melanoplus boulderensis"),], 
       aes(x=Tx90, y=value, colour=elev, group=elev)) + 
  geom_line()+
  #geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(param~1, drop=TRUE, scales="free") +
  xlab("")+ylab("peak of abundance distribution (day)")+
  scale_color_viridis_d()+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")

#low Tx90: early and abundance 


#pdf("Fig2_PhenFits.pdf", height = 10, width = 8)
#grid_arrange_shared_legend(plot.mu, plot.sig, plot.scale, ncol = 1, nrow = 3, position='right')
##plot_grid(plot.mu, plot.scale, nrow=2, rel_heights=c(1,1.4) )
#dev.off()
