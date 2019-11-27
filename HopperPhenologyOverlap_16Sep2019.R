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

#climate data
clim1= read.csv("Clim1Data.csv")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

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

#all species combinations
sp.combs= combn(specs, 2)
#switch focal species
sp.combs= cbind(sp.combs, sp.combs[c(2,1),])

#set up data storage
po= array(data = 0, dim = c(length(sp.combs[1,]), length(years), length(sites),10) ) #add dimension for multiple metrics

#combine species combinations names
sp.combs.n= paste(sp.combs[1,],sp.combs[2,],sep="_")

#add po names
dimnames(po)[[1]] <- sp.combs.n
dimnames(po)[[2]] <- years
dimnames(po)[[3]] <- sites

#----
# plot out densities
#made season warmth a factor
#make cdd_july index
#dat$swarm= as.factor(round(dat$Cdd_siteave,2))
dat$swarm= as.factor(round(dat$Cdd_july_siteave,2))

#adults y=In6
#DIp y= DIp

#plot as proportion of total
dat$elev.lab= paste(dat$elev,"m",sep="")
factor(dat$elev.lab, levels=c("1752m","2195m","2591m","3048m"))
dat$diapause= as.factor(dat$diapause)

#add seasonal timing (calculated below)
dat$timing= timing.mat$timing[match(dat$species, timing.mat$species)]
##drop species without timing
dat= dat[which(!is.na(dat$timing)),]

#order species by timing
dat$species= factor(dat$species, levels=timing.mat$species )

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

## FIGURE 1
#No normalization
pdf("Fig1_distyear_byspec.pdf",height = 12, width = 12)
ggplot(data=dat, aes(x=ordinal, y = DIp, group=spsiteyear, color=Cdd_siteave))+ 
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~elev.lab, scales="free")+
  theme(strip.text.y = element_text(angle = 0))+ theme(legend.position="bottom", legend.key.width=unit(3,"cm"))+
  scale_color_gradientn(colours = c('blue', 'cadetblue', 'orange')) +
  xlab("ordinal date") +ylab("abundance")+ 
   labs(color = "seasonal GDDs")+ theme(strip.text = element_text(face = "italic")) 
dev.off()

#========================================
#Fit distributions
#fit beta: mode, breadth, skewness

#fit gaussian
fitG =
  function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f, method="L-BFGS-B", lower=c(100,20,50), upper=c(300,200,5000) )
  }
#fit= fitG(temp,y,mu=200, sig=1, scale=1)

#data storage
fits= array(data=NA, dim=c(length(specs), length(years), length(sites),3), list(specs,years,sites,c("mu","sig","scale")) )

colors= rainbow(9)

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
      if(p$convergence==0) fits[spec.ind, year.ind, site, 1:3]= p$par
      
      if(year.ind==1 | plot==FALSE) plot(dat.sub2$ordinal, dat.sub2$DIp, type="p", xlim=range(100,300), ylim=range(0,100), col= colors[year.ind], main= paste(sites[site],specs[spec.ind],sep=" ") )
      plot=TRUE
      if(year.ind>1) points(dat.sub2$ordinal, dat.sub2$DIp, type="p", col= colors[year.ind])
      #plot fit
      lines(dat.sub2$ordinal,p$par[3]*dnorm(dat.sub2$ordinal,p$par[1],p$par[2]), col= colors[year.ind])
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

#add seasonal timing (calculated below)
fit.df$timing= timing.mat$timing[match(fit.df$species, timing.mat$species)]
#drop species without timing
fit.df= fit.df[which(!is.na(fit.df$timing)),]

#order species by timing
fit.df$species= factor(fit.df$species, levels=timing.mat$species )

#add elevation label
fit.df$elev.lab= paste(fit.df$elev, "m", sep="")

#cast
rownames(fit.df) <- c()
fit.df1 <- dcast(fit.df, species + year + site + Cdd_siteave ~ param, value.var="value")
fit.df1$spsiteyear= paste(fit.df1$species,fit.df1$year, fit.df1$site, sep="_")

#--------------------
# FIGURE 2

fit.df$param.lab="peak of abundance distribution"
fit.df$param.lab[which(fit.df$param=="sig")]= "breadth of abundance distribution"
fit.df$param.lab[which(fit.df$param=="scale")]= "scale of abundance distribution"
fit.df$param.lab= factor(fit.df$param.lab, levels=c("peak of abundance distribution","breadth of abundance distribution","scale of abundance distribution") )

#set species colors
levels(fit.df$species)
cols= c(brewer.pal(4,"Blues")[2:4],  brewer.pal(5,"Greens")[4:5], brewer.pal(9,"YlOrRd") )
#image(1:14,1,as.matrix(1:14),col=cols)

# ggplot(fit.df[which(fit.df$param=="mu"),], aes(x=cdd_seas, y=value, colour=as.factor(timing), group=species)) +
#   geom_line() + theme_classic() +
#   facet_grid(~elev.lab, drop=TRUE, scales="free") +
#   xlab("")+ylab("peak of abundance distribution")

#DROP shape=period

plot.mu <- ggplot(fit.df[which(fit.df$param=="mu"),], aes(x=cdd_seas, y=value, colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("")+ylab("peak of abundance distribution (day)")+scale_color_manual(values=cols)+ theme(legend.text = element_text(face = "italic")) #+ theme(legend.position="none")

plot.sig <- ggplot(fit.df[which(fit.df$param=="sig"),], aes(x=cdd_seas, y=log(value), colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("")+ylab("log(breadth of abundance distribution) (days)")+scale_color_manual(values=cols)+ theme(legend.text = element_text(face = "italic"))

plot.scale <- ggplot(fit.df[which(fit.df$param=="scale"),], aes(x=cdd_seas, y=log(value), colour=species, group=species)) + geom_point()+geom_smooth(method="lm", se=FALSE)+
  theme_classic() +  #geom_line() +
  facet_grid(~elev.lab, drop=TRUE, scales="free") +
  xlab("seasonal growing degree days (C)")+ylab("log(scale of abundance distribution) (individuals)")+scale_color_manual(values=cols)+ theme(legend.text = element_text(face = "italic"))

pdf("Fig2_PhenFits.pdf", height = 10, width = 8)
grid_arrange_shared_legend(plot.mu, plot.sig, plot.scale, ncol = 1, nrow = 3, position='right')
#plot_grid(plot.mu, plot.scale, nrow=2, rel_heights=c(1,1.4) )
dev.off()

#---------------------
#TABLE STATS
#mu, sig, scale
param.mu= fit.df[which(fit.df$param=="mu"),] #mu scale sig
param.mu= na.omit(param.mu)
param.sig= fit.df[which(fit.df$param=="sig"),] #mu scale sig
param.sig= na.omit(param.sig)
param.scale= fit.df[which(fit.df$param=="scale"),] #mu scale sig
param.scale= na.omit(param.scale)

##omit nymphal diapausers
#param.mu= na.omit(param.mu[which(param.mu$timing!=1),])

#STATS ****
#lme model
mod1=lme(value~ cdd_seas +cdd_seas:elev + cdd_seas:timing + cdd_seas:timing:elev, random=~1|species, data=param.mu )

#TABLE 1
mod.mu=lme(value~ cdd_seas +timing +elev +cdd_seas:timing, random=~1|species, data=param.mu )
mod.sig=lme(value~ cdd_seas+elev, random=~1|species, data=param.sig )
mod.scale=lme(value~ cdd_seas +timing +elev +cdd_seas:timing +elev:timing, random=~1|species, data=param.scale )

summary(mod.sig)
anova(mod.mu)

dredge(mod.scale)

#-------
#plots

#devtools::install_github("hohenstein/remef")
library(remef)

model=lmer(value~ cdd_seas +timing +elev +cdd_seas:timing+(1|species), REML=TRUE, data=param.mu)

y_partial <- remef(model, fix = c("timing", "elev"), ran = "all")
plot(param.mu$cdd_seas, y_partial)
y_partial <- remef(model, fix = c("cdd_seas", "elev"), ran = "all")
plot(param.mu$timing, y_partial)
y_partial <- remef(model, fix = c("cdd_seas", "timing"), ran = "all")
plot(param.mu$elev, y_partial)

#=================================================
#Calculate overlap metrics

#try dropping species
#dat= dat[-which(dat$species=="Melanoplus bivittatus"),]
#dat= dat[-which(dat$species=="Melanoplus sanguinipes"),]

for(year in 1:length(years)){
  for(site in 1:length(sites)){
    
    #sum all individuals by species and year 
    dat.sub= dat[dat$year==years[year],]
    dat.sub= dat.sub[ which(dat.sub$site==sites[site]),]
    #restrict to adults
    #dat.sub= subset(dat.sub, dat.sub$in6>0 )
    
    #change to GDD
    #dat.sub$ordinal= dat.sub$cdd_sumfall
    
    if( nrow(dat.sub)>0 ){ #check data exists
      
      ords= unique(dat.sub$ordinal)
      
      #totals by ordinal date
      dat.agg= dat.sub[,c("species", "in6")] #extract adults
      dat.agg$Species= as.factor(dat.sub$species)
      dat.agg= aggregate(dat.agg$in6, by=list(dat.agg$species), FUN=sum) #sum number adults
      names(dat.agg)= c("Species", "YrTot")
      
      #turn into proportions
      match1= match( dat.sub$species, dat.agg$Species)
      dat.sub$YrTot= dat.agg$YrTot[match1] #Add total number adults per species per year
      dat.sub= dat.sub[which(dat.sub$YrTot>0),]
      dat.sub$AdultProp= dat.sub$in6/dat.sub$YrTot #Add proportion of annual individuals for each ordinal date, normalize species with different abundances 
      #subset to dates with adults
      dat.sub= dat.sub[which(dat.sub$AdultProp>0),]
      
      #find species combinations potentially present
      specs1= unique(dat.sub$species)
      #find species combinations present
      inds=which(sp.combs[1,]%in%specs1 & sp.combs[2,]%in%specs1)
      
      #loop through species combinations
      for(ind.k in 1:length(inds) ){
        sp1= sp.combs[1,inds[ind.k] ]
        sp2= sp.combs[2,inds[ind.k] ]
        #subset to combinations
        dat.sub2= dat.sub[which(dat.sub$species==sp1 | dat.sub$species==sp2),]
        
        #   #Pianka (1974) metric, http://www.jstor.org/stable/pdf/4217327.pdf, 
        
        #new metric?: sum across dates: (max proportion-overlap)/total, On the Analysis of Phenological Overlap
        #   #po[inds[c], year, site]= po[inds[c], year, site]+ (max(prop1,prop2)-abs(prop1-prop2))/(prop1+prop2)
       
        #DI overlap
       
        #check each species has at least 4 dates
        if( length(which(dat.sub2$species==sp1))>3 & length(which(dat.sub2$species==sp2))>3 ){
        
        #loess curve for each species
        d.sp1= dat.sub2[which(dat.sub2$species==sp1),]
        d.sp2= dat.sub2[which(dat.sub2$species==sp2),]
        
        # create loess object and prediction function
        l1 <- loess(DIp~ordinal, d.sp1)
        l2 <- loess(DIp~ordinal, d.sp2)
        
        #set negative values and NA to zero
        f1 <- function(x) {f=predict(l1,newdata=x); f[f<0]=0;f[is.na(f)]=0; return(f)}
        f2 <- function(x) {f=predict(l2,newdata=x); f[f<0]=0;f[is.na(f)]=0; return(f)}
        
        #find range
        low= min(dat.sub2$ordinal)
        up= max(dat.sub2$ordinal)
       
        # perform integration
        a1=integrate(f1,low,up)$value
        a2=integrate(f2,low,up)$value
     
        #area overlap
        #normalize to 1
        sp1.n= f1(low:up)/a1
        sp2.n= f2(low:up)/a2
        #pianka's index
        po[inds[ind.k], year, site,1]= sum(sp1.n*sp2.n)/sqrt(sum(sp1.n^2)*sum(sp2.n^2))
        
        #overlapping area normalized to total area
        #normalize to 1
        sp1.n= f1(low:up)
        sp2.n= f2(low:up)
        dinds= which(sp1.n>0 & sp2.n>0)
        if(length(dinds>0)) po[inds[ind.k], year, site,2]= sum(sapply(dinds, FUN=function(x) min(sp1.n[x],sp2.n[x])) )/(a1+a2)
        #overlapping area normalized to total area of focal species
        if(length(dinds>0)) po[inds[ind.k], year, site,3]= sum(sapply(dinds, FUN=function(x) min(sp1.n[x],sp2.n[x])) )/(a1)
        #not normalized
        if(length(dinds>0)) po[inds[ind.k], year, site,4]= sum(sapply(dinds, FUN=function(x) min(sp1.n[x],sp2.n[x])) )
        
        #area of focal species
        po[inds[ind.k], year, site,5]= a1
        
        #time difference metrics
        #days focal
        po[inds[ind.k], year, site,6]= length(which(sp1.n>0))
        #days overlap
        po[inds[ind.k], year, site,7]= length(dinds)
        #days overlap normalized
        po[inds[ind.k], year, site,8]= length(dinds)/length(which(sp1.n>0))
        #days difference in peak
        po[inds[ind.k], year, site,9]= abs(which.max(sp1.n)-which.max(sp2.n))
        #days difference in 20% of peak
        po[inds[ind.k], year, site,10]= abs(which.max(sp1.n>max(sp1.n)*0.2) - which.max(sp2.n>max(sp2.n)*0.2) )
        
        } #end check each species has at least 4 dates
      } #end loop species combinations
      
    } #end check data
  } #end loop year
} #end loop sites

#replace 0 values with NA
po[po==0]=NA

#-----------------------------------
#PLOTS
# PO: sp comb x years x sites

po1= melt(po, varnames=c("sp","year","site","metric"))

#drop NAs
po1= po1[which(!is.na(po1$value)),]

#add time period
po1$period="initial"
po1[which(po1$year>1960),"period"]<-"resurvey"

#ADD cdd and temp
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
po1$siteyear= paste(po1$site, po1$year, sep="")

#select cdd metric
clim1$Cdd= clim1$Cdd_seas

po1$Tmean=NA
po1$cdd=NA
po1$cdd_july=NA
match1= match(po1$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean[matched]<- clim1$Mean[match1[matched]]  
po1$cdd[matched]<- clim1$Cdd[match1[matched]]
po1$cdd_july[matched]<- clim1$Cdd_july[match1[matched]]

#---
#Change years for plotting ### FIX
po1[which(po1$year==1958),"year"]= 1958+40
po1[which(po1$year==1959),"year"]= 1959+40
po1[which(po1$year==1960),"year"]= 1960+40

#separate species names for plotting
po1$sp= as.character(po1$sp)

sps= colsplit(po1$sp, "_", names = c('sp1', 'sp2'))
po1$sp1<- sps$sp1
po1$sp2<- sps$sp2

#order by early and late species
#ave phen
#hop.el = hop1 %>% group_by(species,site) %>% arrange(species,site) %>% mutate(phen = mean(ordinal))
hop.agg= aggregate(hop1, list(hop1$species),FUN=mean) #fix for hop1$site,
hop.agg= hop.agg[order(hop.agg$ordinal),]

#make species factor for plotting
po1$sp1= factor(po1$sp1, levels=specs)
po1$sp2= factor(po1$sp2, levels=specs)

rownames(timing.mat)<- timing.mat$species
timing.mat$fits1= as.numeric(as.character(timing.mat$fits1))

#add seasonal timing (calculated below)
po1$sp1_timing= timing.mat$timing[match(po1$sp1, rownames(timing.mat))]
po1$sp2_timing= timing.mat$timing[match(po1$sp2, rownames(timing.mat))]
po1$sp1_timing_doy= round(timing.mat$fits1[match(po1$sp1, rownames(timing.mat))])
po1$sp2_timing_doy= round(timing.mat$fits1[match(po1$sp2, rownames(timing.mat))])

##drop species without timing
po1= po1[which(!is.na(po1$sp1_timing)),]
po1= po1[which(!is.na(po1$sp2_timing)),]

#PLOT
#add elevation
po1$elevation= as.factor(elevs[match(po1$site, sites)])

#*********************************************
#Plot overlap

#assess significance
#calculate significant regressions
#apply through combinations of species and elevations
po1$elevspec= paste(po1$elevation, po1$sp, sep="")
elevspec= matrix(unique(po1$elevspec))

#add average year CDD
fit.df$siteyear= paste(fit.df$site, fit.df$year, sep="")
po1$Cdd_siteave = fit.df[match(po1$siteyear, fit.df$siteyear),"Cdd_siteave"]

#---
#EXTRACT COEFFS
metric.ind=3

#2: overlap area normalized to 1
#7: days overlap
#9: days difference in peak
metrics=c(2,7,9)
po2= po1[which(po1$metric==metrics[metric.ind]),]
x1="cdd" # "Cdd_siteave", "cdd", "cdd_july"
#x1="year"
#x1="Tmean"

po2$elevation= as.character(po2$elevation)

p.gdd= apply(elevspec,1, FUN=function(x) tryCatch( summary(lm(value~get(x1), data=po2[which(po2$elevation==substr(x,1,4)& po2$sp==substr(x,5,nchar(x)) ),] ))$coefficients[2,], error=function(err) rep(NA,4)))

#combine
p.mat=as.data.frame(cbind(elevspec,t(p.gdd) ))
colnames(p.mat)[2:5]= c("Estimate","Std. Error","t value","Pr(>|t|)") 

#make numeric
p.mat$Estimate= as.numeric(as.character(p.mat$Estimate))
p.mat$`Std. Error`= as.numeric(as.character(p.mat$`Std. Error`))
#CIs
p.mat$upCI= p.mat$Estimate + 1.96*p.mat$`Std. Error`
p.mat$lowCI= p.mat$Estimate - 1.96*p.mat$`Std. Error`

#add columns for significance
p.mat$sig.gdd<-"nonsignificant"
p.mat$sig.gdd[which(p.gdd[4,]<0.05)]="significant"

#add back to matrix
match1= match(po2$elevspec, elevspec)
po2$sig.gdd= factor(p.mat[match1,"sig.gdd"], levels=c("significant","nonsignificant"))

#---

#PLOTS

#drop if less than 5 years 
po2.counts= ddply(po2, c("site", "sp","elevation","elevspec"), summarise, count=length(unique(year)) )
po2.keep= po2.counts[which(po2.counts$count>3),]
po2= po2[which(po2$elevspec %in% unique(po2.keep$elevspec)),]

#order species by timing
po2$sp1= factor(po2$sp1, levels=row.names(timing.mat) )

#labels
#add elevation label
po2$elev.lab= paste(po2$elevation, "m", sep="")
metric.lab= c("overlapping area (proportion)", "days of overlap", "days difference in peak abundance" )

#plot by elevation
if(metric.ind==1) pdf("Fig3S_PhenOverlap_byGDD_overlaparea.pdf", height = 10, width = 10)
if(metric.ind==2) pdf("Fig3S_PhenOverlap_byGDD_daysoverlap.pdf", height = 10, width = 10)
if(metric.ind==3) pdf("Fig3S_PhenOverlap_byGDD_daysdiffpeak.pdf", height = 10, width = 10)
ggplot(data=po2, aes_string(x=x1, y = "value", group="sp2", color="sp2_timing_doy"))+geom_point(aes(shape=period, fill=period), size=2)+
  facet_grid(sp1~elev.lab, drop=TRUE, scales="free_x")+theme_bw()+geom_smooth(method="lm", se=FALSE)+
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_color_gradientn(colours = c('blue', 'cadetblue', 'orange')  ) +
  xlab("seasonal growing degree days") +ylab(metric.lab[metric.ind])+ 
  theme(strip.text = element_text(face = "italic")) +
  labs(color = "seasonal timing \n (doy)")
#, aes(linetype=sig.gdd)
dev.off()

#+ylab("Phenological overlap (days)")+xlab("Growing degree days")

#------
#STATS ACROSS SPECIES #*****
po2$elevation= as.numeric(po2$elevation)

mod1=lme(value~ cdd*sp1_timing*sp2_timing*elevation, random=~1|sp1/sp2, data=po2)

mod1=lme(value~ cdd + cdd:elevation + cdd:sp1_timing +cdd:sp2_timing+ cdd:sp1_timing:elevation +cdd:sp2_timing:elevation+cdd:sp1_timing:sp2_timing, random=~1|sp1/sp2, data=po2)

mod1=lme(value~ cdd +elevation +sp1_timing +sp2_timing +cdd:elevation + cdd:sp1_timing +cdd:sp2_timing +elevation:sp1_timing  +elevation:sp2_timing +sp1_timing:sp2_timing +cdd:sp1_timing:sp2_timing +elevation:sp1_timing:sp2_timing, random=~1|sp1/sp2, data=po2)

#dredge(mod1)
summary(mod1)
anova(mod1)

#-------------
#STATS

#MAKE TIMING MAT
#C1 timing
# fits1= fits[,,3,1]
# fits1= sort(rowMeans(fits1, na.rm=TRUE))
# 
# #Calcualte average timing
# #early vs late species
# fits1= fits[,,2,1]
# fits1= sort(rowMeans(fits1, na.rm=TRUE))
# timing.mat= as.data.frame( cbind(colnames(fits1), fits1) )
# timing.mat$timing= c(rep(1,3),rep(2,2),rep(3,8) )
# #Drop Trimerotropis cincta due to limited data
# #timing.mat= timing.mat[1:13,]
# #ADD M. fasciatus 
# timing.mat= rbind(timing.mat, c(NA, 3))
# rownames(timing.mat)[nrow(timing.mat)]<- "Melanoplus fasciatus"
# 
# #write out
# setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
# #write.csv(timing.mat, "timingmat_2019.csv" )
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")
#-------------------

#pick species
po2.sp= po2[which(po2$sp1==rownames(timing.mat)[7]),]
po2.sp$timing= timing.mat$timing[match(po2.sp$sp2, rownames(timing.mat))]
#drop species without timing
po2.sp= na.omit(po2.sp[which(!is.na(po2.sp)),] )

#make elevation numeric
po2.sp$elevation= as.numeric(as.character(po2.sp$elevation))

mod1=lme(value~ cdd + elevation + timing, random=~1|sp2, data=po2.sp)
print(po2.sp$sp1[1] )
summary(mod1)

#===================================
# STATS GROUP BY SPECIES
po2= po1[which(po1$metric==4),] 

mod1= lm(value~ cdd*elevation+sp2+sp1, data=po2)
mod1= lm(value~ cdd*as.numeric(as.character(elevation))+sp, data=po2)

#---------------
#save overlap output in table

overlaps.sp= matrix(NA, nrow=length(specs), ncol=3 )
rownames(overlaps.sp)= specs

#by focal species

for(sp.k in 1:length(specs)){
po3= subset(po2, po2$sp1==specs[sp.k] ) #subset(po2, po2$sp1==specs[sp.k] | po2$sp2==specs[sp.k])
# #switch species order
# p.temp= po3$sp1
# inds= which(po3$sp2==specs[sp.k])
# if(length(inds)>0){
#  po3$sp1[inds]=specs[sp.k] 
#  po3$sp2[inds]=po3$sp1[inds] 
# }

if(nrow(po3)>0) {

po3$elevation= as.numeric(as.character(po3$elevation))

#---------------
#early vs late
po3$timing<-"early"
po3$timing[which(po3$sp2 %in% specs[c(1,4)])]<-"late"

overlaps.sp[sp.k,]= tryCatch( summary(lme(value~ cdd_july+cdd_july:elevation+cdd_july:timing,random=~1|sp2 , data=po3))$tTable[2:4,5], error=function(err) rep(NA,3))
#mod1= lme(value~ cdd_july+elevation+timing+cdd_july:elevation+cdd_july:timing,random=~1|sp2 , data=po3)
} #end data check
} #end species loop

#======================
#Average overlap across community

#drop species combinations with less that 4 years
po1.counts= ddply(po1, c("site", "sp","elevation","elevspec","metric"), summarise, count=length(unique(year)) )
po1.keep= po1.counts[which(po1.counts$count>3),]
po1.keep$elevspecmetric= paste(po1.keep$elevspec, po1.keep$metric , sep="")
po1$elevspecmetric= paste(po1$elevspec, po1$metric , sep="")

po1= po1[which(po1$elevspecmetric %in% unique(po1.keep$elevspecmetric)),]

#-----------
#AVERAGE
#Average across seasonal timing
po1$sp1_sp2= paste(po1$sp1_timing, po1$sp2_timing, sep="_")

#by timing groups
po.tim= ddply(po1, c("site", "year","metric","period","siteyear","sp1_sp2","sp1_timing", "sp2_timing"), summarise,
               value = mean(value, na.rm=TRUE), Tmean= Tmean[1], cdd=cdd[1],cdd_july=cdd_july[1],elevation=elevation[1], Cdd_siteave= Cdd_siteave[1])
po.tim$elev.lab= paste(po.tim$elevation,"m",sep="")
#symmetric, so drop repeated combinations
po.tim<- po.tim[-which(po.tim$sp1_sp2 %in% c("2_1","3_1","3_2")),]

#by timing
po.st= ddply(po1, c("site", "year","metric","period","siteyear","sp1","sp1_timing","sp2_timing"), summarise,
              value = mean(value, na.rm=TRUE), Tmean= Tmean[1], cdd=cdd[1],cdd_july=cdd_july[1],elevation=elevation[1], Cdd_siteave= Cdd_siteave[1])
po.st$elev.lab= paste(po.st$elevation,"m",sep="")

#Average across comminuty
po.comm= ddply(po1, c("site", "year","metric","period","siteyear"), summarise,
               value = mean(value, na.rm=TRUE), Tmean= Tmean[1], cdd=cdd[1],cdd_july=cdd_july[1],elevation=elevation[1], Cdd_siteave= Cdd_siteave[1])
po.comm$elev.lab= paste(po.comm$elevation,"m",sep="")

#------------
#FIGURE : plot by seasonal timing
#How to color background: https://stackoverflow.com/questions/6750664/how-to-change-the-format-of-an-individual-facet-wrap-panel

#timing * species

#order species by timing
po.st$sp1= factor(po.st$sp1, levels=row.names(timing.mat) )
po.st$st.lab= NA
po.st$st.lab[po.st$sp2_timing==1]<-"nymphal diapausers"
po.st$st.lab[po.st$sp2_timing==2]<-"early season"
po.st$st.lab[po.st$sp2_timing==3]<-"late season"
po.st$st.lab= factor(po.st$st.lab, levels=c("nymphal diapausers","early season","late season") )

pdf("FigS1_CommOverlap_bysptiming_m2.pdf", height = 10, width = 10)
ggplot(data=po.st[which(po.st$metric==2),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+facet_grid(sp1~st.lab, drop=TRUE, scales="free_x") +theme(strip.text.y = element_text(angle = 0))+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"), name = "Elevation (m)") +theme(legend.position="bottom")+
  ylab("overlapping area (proportion)")+xlab("seasonal growing degree days (C)")+ 
  theme(strip.text = element_text(face = "italic")) 
dev.off()

pdf("FigS2_CommOverlap_bysptiming_m7.pdf", height = 10, width = 10)
ggplot(data=po.st[which(po.st$metric==7),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+facet_grid(sp1~st.lab, drop=TRUE, scales="free_x") +theme(strip.text.y = element_text(angle = 0))+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"), name = "Elevation (m)") +theme(legend.position="bottom")+
  ylab("days of overlap")+xlab("seasonal growing degree days (C)")+ 
  theme(strip.text = element_text(face = "italic")) 
dev.off()

pdf("FigS3_CommOverlap_bysptiming_m9.pdf", height = 10, width = 10)
ggplot(data=po.st[which(po.st$metric==9),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+facet_grid(sp1~st.lab, drop=TRUE, scales="free_x") +theme(strip.text.y = element_text(angle = 0))+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"), name = "Elevation (m)") +theme(legend.position="bottom")+
  ylab("days difference in peak abundance")+xlab("seasonal growing degree days (C)")+ 
  theme(strip.text = element_text(face = "italic")) 
dev.off()

#-----
po.tim$fst.lab= NA
po.tim$fst.lab[po.tim$sp1_timing==1]<-"focal: nymphal diapausers"
po.tim$fst.lab[po.tim$sp1_timing==2]<-"focal: early season"
po.tim$fst.lab[po.tim$sp1_timing==3]<-"focal: late season"
po.tim$st.lab= NA
po.tim$st.lab[po.tim$sp2_timing==1]<-"compared: nymphal diapausers"
po.tim$st.lab[po.tim$sp2_timing==2]<-"compared: early season"
po.tim$st.lab[po.tim$sp2_timing==3]<-"compared: late season"

po.tim$fst.lab= factor(po.tim$fst.lab, levels=c("focal: nymphal diapausers","focal: early season","focal: late season") )
po.tim$st.lab= factor(po.tim$st.lab, levels=c("compared: nymphal diapausers","compared: early season","compared: late season") )

#FIG 3
pdf("Fig3_CommOverlap_bytiming.pdf", height = 10, width = 10)
ggplot(data=po.tim[which(po.tim$metric==2),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+facet_wrap(fst.lab~st.lab, drop=TRUE, scales="free_x")+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3")) +
  theme(legend.position="bottom")+ylab(metric.lab[1])+xlab("seasonal growing degree days (C)")+ 
  guides(color=guide_legend(title="elevation (m)" ))
dev.off()

#TABLE 2 STATS
#STATS ACROSS SPECIES #***** #metrics 2,7,9
po.tim2= na.omit(po.tim[which(po.tim$metric==9),])
po.tim2$elevation= as.numeric(as.character(po.tim2$elevation))

mod1=lm(value~ cdd+sp1_timing+sp2_timing+elevation+cdd:sp1_timing+cdd:sp2_timing+cdd:elevation+sp1_timing:sp2_timing, data=po.tim2)

#dredge(mod1)
summary(mod1)
anova(mod1)

#--------------------------
#partial residual plot
library(faraway)

my.prplot= function (g, i, xlabs)  #beautify plot
{
  xl <- xlabs[i]
  yl <- paste("")
  x <- model.matrix(g)[, i + 1]
  plot(x, g$coeff[i + 1] * x + g$res, xlab = xl, ylab = yl, col=rgb(0.2,0.2, 0.2, 0.5),pch=16)        #@I played with different plot characters to get different plots. Not sold
  abline(0, g$coeff[i + 1])
  invisible()
}
#--------------

xlabs= c("GDDs","focal species timing", "compared species timing", "elevation (m)", "GDDs:focal species timing", "GDDs:compared species timing", "GDDs:elevation (m)", "focal:compared species timing")

pdf("FigS4_OverlapPR.pdf",height = 10, width = 10)

par(mfrow=c(3,3), cex=1.1, lwd=1, mar=c(3,2, 1, 1), mgp=c(1.3, 0.5, 0), oma=c(0,2,0,0), bty="l", cex.lab=1.2)
my.prplot(mod1, 1, xlabs)
my.prplot(mod1, 2, xlabs)
my.prplot(mod1, 3, xlabs)
my.prplot(mod1, 4, xlabs)
my.prplot(mod1, 5, xlabs)
my.prplot(mod1, 6, xlabs)
my.prplot(mod1, 7, xlabs)
my.prplot(mod1, 8, xlabs)

mtext("Partial residual of overlapping area (proportion)", side=2, line = -0.5, cex=1.3, outer=TRUE)

dev.off()

#------------
#FIGURE 4: plot accross community
#by year
#ggplot(data=po.comm[which(po.comm$metric==3),], aes(x=year, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+theme_bw()+geom_smooth(method="lm", se=FALSE)

#all metrics
plot.m1= ggplot(data=po.comm[which(po.comm$metric==1),], aes(x=cdd, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+labs(title="Pianka index")+ theme(legend.position="none")
plot.m3= ggplot(data=po.comm[which(po.comm$metric==3),], aes(x=cdd, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+labs(title="overlap area, norm to area focal")+ theme(legend.position="none")
plot.m8= ggplot(data=po.comm[which(po.comm$metric==8),], aes(x=cdd, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)+labs(title="days overlap, norm")+ theme(legend.position="none")

#focus on 2,7,9?
plot.m2= ggplot(data=po.comm[which(po.comm$metric==2),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE) +theme(legend.position="none") +ylab(metric.lab[1])+xlab("")+
 scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))+ 
  guides(color=guide_legend(title="elevation (m)" ))

plot.m7= ggplot(data=po.comm[which(po.comm$metric==7),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE) +theme(legend.position="none") +ylab(metric.lab[2])+xlab("")+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))+ 
  guides(color=guide_legend(title="elevation (m)" ))

plot.m9= ggplot(data=po.comm[which(po.comm$metric==9),], aes(x=cdd, y = value, color=elev.lab))+geom_point(aes(), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE) +theme(legend.position="bottom") +ylab(metric.lab[3])+xlab("seasonal growing degree days (C)")+
  scale_color_manual(breaks = c("1752m", "2195m", "2591m","3048m"),
                     values=c("darkorange3", "darkorange", "cornflowerblue","blue3"))+ 
  guides(color=guide_legend(title="elevation (m)" ))

#FIGURE 4
pdf("Fig4_CommOverlap.pdf", height = 10, width = 6)
plot_grid(plot.m2, plot.m7, plot.m9, nrow=3, labels=c("a.","b.","c."), label_x=0.1)
dev.off()

#----
#STATS *****
po.comm$elevation= as.numeric(as.character(po.comm$elevation))

mod1= lm(value~ cdd+ cdd:elevation, data=po.comm[which(po.comm$metric==2),])
mod1= lm(value~ cdd+ cdd:elevation, data=po.comm[which(po.comm$metric==7),])
mod1= lm(value~ cdd+ cdd:elevation, data=po.comm[which(po.comm$metric==9),])
summary(mod1)
anova(mod1)

prplot(mod1, 1)
prplot(mod1, 2)

#single elevation
mod1= lm(value~ cdd, data=po.comm[which(po.comm$metric==9 & po.comm$elevation==3048),])

#library(faraway)
prplot(mod1,1)
prplot(mod1,2)


