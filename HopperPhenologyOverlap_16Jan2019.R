library(reshape)
library(reshape2)
library(tidyr)

#LOAD SPECIES DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/data/")
spec.diapause= read.csv("SpeciesDiapause.csv")
dat.all= read.csv("HopperClimateData.csv")
#data currently from phenology file PhenFigs_20May2018.R

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

#------------------------------------
#compare egg diapausers vs non
dat.all$diapause="egg"
dat.all$species= as.character(dat.all$species)
dat.all$diapause[which(dat.all$species %in% c("Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes"))]="nymph"

#fix names
dat.all$species[which(dat.all$species=="Melanoplus bivatattus")]="Melanoplus bivittatus"

#focal species
specs= c("Aeropedellus clavatus","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes","Arphia conspersa","Eritettix simplex","Pardalaphoa apiculata","Xanthippus corallipes") #"Chloealtis abdominalis",
#subset to focal species
##dat= subset(dat.all, dat.all$species %in% specs)
#keep all species
dat=dat.all

years= unique(na.omit(hop$year))
sites= unique(hop$site)

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

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

pdf("Fig_dipnorm.pdf",height = 11, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIpNorm, group=species, color=species, lty=diapause))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(year~elev.lab, scales="free")+ylim(0,0.2)+
  ylab("Proportional abundance")+xlab("Day of year")+ theme(legend.position="bottom")+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
dev.off()

pdf("Fig_dipnormyear.pdf",height = 10, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIpNorm, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(year~elev.lab, scales="free")+ylim(0,0.2)+ylab("Proportional abundance")+xlab("Day of year")
dev.off()

#---
#by species
pdf("Fig_distnormyear_byspec.pdf",height = 20, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIpNorm, group=spsiteyear, color=Cdd_siteave, lty=diapause))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~site, scales="free")+ylim(0,0.2)+
  theme(strip.text.y = element_text(angle = 0))
dev.off()

#No normalization
pdf("Fig_distyear_byspec.pdf",height = 20, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIp, group=spsiteyear, color=Cdd_siteave, lty=diapause))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~site, scales="free")+
  theme(strip.text.y = element_text(angle = 0))
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
    
    optim(c(mu,sig,scale),f)
  }
#fit= fitG(temp,y,mu=200, sig=1, scale=1)

#data storage
fits= array(data=NA, dim=c(length(specs), length(years), length(sites),3), list(specs,years,sites,c("mu","sig","scale")) )

colors= rainbow(9)

pdf("Fits.pdf", height = 10, width = 10)
par(mfrow=c(5,5))

for(year in 1:length(years)){
  for(site in 1:length(sites)){
    
    dat.sub= dat[which(dat$year==years[year] & dat$site==sites[site] ),]

    spec.inds= which(specs %in% dat.sub$species )
    
    for(spec.ind in spec.inds){
      dat.sub2= dat.sub[which(dat.sub$species==specs[spec.ind]),]
  
      p= fitG(dat.sub2$ordinal, dat.sub2$DIp, mu=150, sig=20, scale=100 )
      if(p$convergence==0) fits[spec.ind, year, site, 1:3]= p$par
      
      if(spec.ind==1) plot(dat.sub2$ordinal, dat.sub2$DIp, type="p", xlim=range(100,300), ylim=range(0,100), col= colors[spec.ind])
      if(spec.ind>1) points(dat.sub2$ordinal, dat.sub2$DIp, type="p", col= colors[spec.ind])
      #plot fit
      lines(dat.sub2$ordinal,p$par[3]*dnorm(dat.sub2$ordinal,p$par[1],p$par[2]), col= colors[spec.ind])
      
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

#plot
#mu, sig, scale
pdf("PhenFits.pdf", height = 20, width = 10)
ggplot(data=fit.df[which(fit.df$param=="mu"),], aes(x=Cdd_siteave, y = value, group=elev, color=as.factor(elev)))+
  geom_point(aes(shape=period, fill=period))+facet_grid(species~., scales="free")+geom_smooth(method="lm", se=FALSE)+geom_line()
dev.off()

#stats
#mu
param.mu= fit.df[which(fit.df$param=="mu"),]
param.mu= param.mu[which(param.mu$value>0),]
param.mu= na.omit(param.mu)

mod1= lm(value~ species+elev+year, data=param.mu)

#=================================================
#Calculate overlap metrics

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
#clim1$Cdd=  clim1$Cdd_sum
#clim1$Cdd=  clim1$Cdd_june
#clim1$Cdd=  clim1$Cdd_july
#clim1$Cdd=  clim1$Cdd_aug
#clim1$Cdd=  clim1$Cdd_early
#clim1$Cdd=  clim1$Cdd_mid
clim1$Cdd= clim1$Cdd_seas

po1$Tmean=NA
po1$cdd=NA
po1$cdd_july=NA
match1= match(po1$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean[matched]<- clim1$Mean[match1[matched]]  
po1$cdd[matched]<- clim1$Cdd[match1[matched]]
po1$cdd_july[matched]<- clim1$Cdd_july[match1[matched]]

# ##match back to dat
# datm= dat[duplicated(dat$siteyear)==FALSE,]
# match1= match(po1$siteyear, datm$siteyear)
# matched= which(!is.na(match1))
# po1$Tmean[matched]<- datm$Tmean[match1[matched]]  
# po1$cdd[matched]<- datm$cdd_seas[match1[matched]] 

#---
#Change years for plotting ### FIX
po1[which(po1$year==1958),"year"]= 1958+40
po1[which(po1$year==1959),"year"]= 1959+40
po1[which(po1$year==1960),"year"]= 1960+40

#separate species names for plotting
po1$sp= as.character(po1$sp)

sps= colsplit(po1$sp, pattern = "_", names = c('sp1', 'sp2'))
po1$sp1<- sps$sp1
po1$sp2<- sps$sp2

#order by early and late species
#ave phen
#hop.el = hop1 %>% group_by(species,site) %>% arrange(species,site) %>% mutate(phen = mean(ordinal))
hop.agg= aggregate(hop1, list(hop1$species),FUN=mean) #fix for hop1$site,
hop.agg= hop.agg[order(hop.agg$ordinal),]

#make species factor for plotting
#po1$sp1= factor(po1$sp1, levels=hop.agg$Group.1)
#po1$sp2= factor(po1$sp2, levels=hop.agg$Group.1)
po1$sp1= factor(po1$sp1, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Camnula pellucida","Melanoplus sanguinipes", "Melanoplus dawsoni", "Chloealtis abdominalis"))
po1$sp2= factor(po1$sp2, levels=c("Aeropedellus clavatus","Melanoplus boulderensis","Camnula pellucida","Melanoplus sanguinipes", "Melanoplus dawsoni", "Chloealtis abdominalis"))

#PLOT
#add elevation
po1$elevation= as.factor(elevs[match(po1$site, sites)])

#drop Chataqua for now 
#po1= po1[which(po1$site!="CHA"),]

#focus on sites B1 and C1 for now
#po1= po1[which(po1$site %in% c("B1", "C1")) ,] #[which(po1$metric==2),] #

#*********************************************
#Plot overlap

#assess significance
#calculate significant regressions
#apply through combinations of species and elevations
po1$elevspec= paste(po1$elevation, po1$sp, sep="")
elevspec= matrix(unique(po1$elevspec))

#---
#EXTRACT COEFFS
po2= po1[which(po1$metric==7),] #3,7
x1="cdd_july"
#x1="year"
#x1="Tmean"

p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(value~get(x1), data=po2[which(po2$elevation==substr(x,1,4)& po2$sp==substr(x,5,nchar(x)) ),] ))$coefficients[2,])

#combine
p.mat=as.data.frame(cbind(elevspec,t(p.gdd) ))
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

pdf("PhenOverlap_byGDD_days.pdf", height = 8, width = 10)
ggplot(data=po2, aes_string(x=x1, y = "value", color="elevation"))+geom_point(aes(shape=period, fill=period), size=2)+
  facet_grid(sp1~sp2, drop=TRUE)+theme_bw()+geom_smooth(method="lm", se=FALSE, aes(linetype=sig.gdd))+ 
  scale_color_manual(values=c("darkorange", "blue","darkgreen","purple"))#+xlim(200,750)
dev.off()

#+ylab("Phenological overlap (days)")+xlab("Growing degree days")

#**********************************************

# #overlap by GDD by site
# pdf("PhenOverlap_byGDD_site.pdf", height = 8, width = 10)
# ggplot(data=po1[which(po1$metric==1),], aes(x=cdd, y = value, color=sp))+geom_point(aes(shape=period, fill=period), size=2)+
#   facet_wrap(~site, drop=TRUE, scales="free")+theme_bw()+geom_smooth(method="lm", se=FALSE)+ylab("Phenological overlap")+xlab("Growing degree days")+ 
#   scale_color_manual(values=rainbow(20))+ylim(0,1)
# dev.off()

#--------------
#Average overlap across community

po.comm= ddply(po1, c("site", "year","metric","period","siteyear"), summarise,
               value = mean(value, na.rm=TRUE), Tmean= Tmean[1], cdd=cdd[1],cdd_july=cdd_july[1],elevation=elevation[1])
                 
ggplot(data=po.comm[which(po.comm$metric==1),], aes(x=year, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+
  theme_bw()+geom_smooth(method="lm", se=FALSE)

#===================================
# STATS GROUP BY SPECIES
po2= po1[which(po1$metric==4),] 

mod1= lm(value~ cdd*elevation+sp2+sp1, data=po2)
mod1= lm(value~ cdd*as.numeric(as.character(elevation))+sp, data=po2)

#---------------
#by focal species
sp.k=5
po3= subset(po2, po2$sp1==specs[sp.k] ) #subset(po2, po2$sp1==specs[sp.k] | po2$sp2==specs[sp.k])
# #switch species order
# p.temp= po3$sp1
# inds= which(po3$sp2==specs[sp.k])
# if(length(inds)>0){
#  po3$sp1[inds]=specs[sp.k] 
#  po3$sp2[inds]=po3$sp1[inds] 
# }

po3$elevation= as.numeric(as.character(po3$elevation))
mod1= lm(value~ cdd_july+cdd_july:elevation+cdd_july:sp2, data=po3)
specs[sp.k]
summary(mod1)

#---------------
#early vs late
require(nlme)
require(lme4)

po3$timing="early"
po3$timing[which(po3$sp2 %in% specs[c(1,4)])]<-"late"

#mod1= lme(value~ cdd_july+cdd_july:elevation,random=~1|sp2 , data=po3)
mod1= lme(value~ cdd_july+cdd_july:elevation+cdd_july:timing,random=~1|sp2 , data=po3)
summary(mod1)
#--------------------------------


