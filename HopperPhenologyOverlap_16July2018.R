library(reshape)
library(reshape2)
library(tidyr)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/GrasshopperPhenSynch/figures/")

#data currently from phenology file PhenFigs_20May2018.R
#------------------------------------
dat.all= dat

#focal species
specs= c("Aeropedellus clavatus","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes") #"Chloealtis abdominalis",

years= unique(na.omit(hop$year))
sites= unique(hop$site)

#elevs and sites
sites=c("CHA", "A1", "B1", "C1", "D1")
elevs= c(1752, 2195, 2591, 3048, 3739)

#subset ot focal species
dat= subset(dat, dat$species %in% specs)

#Calculate development index
dat$DIp=0
inds=which(dat$total>0)  
dat$DIp[inds]= (dat$in1[inds] +dat$in2[inds]*2 +dat$in3[inds]*3 +dat$in4[inds]*4 +dat$in5[inds]*5 +dat$in6[inds]*6)/6

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

#set up data storage
po= array(data = 0, dim = c(length(sp.combs[1,]), length(years), length(sites),5) ) #add dimension for multiple metrics

#combine species combinations names
sp.combs.n= paste(sp.combs[1,],sp.combs[2,],sep="_")

#add po names
dimnames(po)[[1]] <- sp.combs.n
dimnames(po)[[2]] <- years
dimnames(po)[[3]] <- sites

#----
# plot out densities
#made season warmth a factor
dat$swarm= as.factor(round(dat$Cdd_siteave,2))

#adults
pdf("Fig_dist.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = in6, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(swarm~site, scales="free")
dev.off()

pdf("Fig_distyear.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = in6, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(year~site, scales="free")
dev.off()

#DIp
pdf("Fig_dip.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIp, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(swarm~site, scales="free")
dev.off()

pdf("Fig_dipyear.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIp, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(year~site, scales="free")
dev.off()

#----
#plot as proportion of total

pdf("Fig_dipnorm.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIpNorm, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(swarm~site, scales="free")+ylim(0,0.2)
dev.off()

pdf("Fig_dipnormyear.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = DIpNorm, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(year~site, scales="free")+ylim(0,0.2)
dev.off()

pdf("Fig_distnorm.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = AdultNorm, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(swarm~site, scales="free")+ylim(0,0.2)
dev.off()

pdf("Fig_distnormyear.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = AdultNorm, group=species, color=species))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(year~site, scales="free")+ylim(0,0.2)
dev.off()

#---
#by species
pdf("Fig_distnormyear_byspec.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = AdultNorm, group=spsiteyear, color=Cdd_siteave))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~site, scales="free")+ylim(0,0.2)
dev.off()

pdf("Fig_distyear_byspec.pdf",height = 12, width = 10)
ggplot(data=dat, aes(x=ordinal, y = in6, group=spsiteyear, color=Cdd_siteave))+
  geom_smooth(method="loess", se=FALSE)+geom_point()+facet_grid(species~site, scales="free")
dev.off()

#----------------------------
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
        
        #time difference metrics
        #days overlap
        po[inds[ind.k], year, site,3]= length(dinds)
        #days difference in peak
        po[inds[ind.k], year, site,4]= abs(which.max(sp1.n)-which.max(sp2.n))
        #days difference in 20% of peak
        po[inds[ind.k], year, site,5]= abs(which.max(sp1.n>max(sp1.n)*0.2) - which.max(sp2.n>max(sp2.n)*0.2) )
        
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
match1= match(po1$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean[matched]<- clim1$Mean[match1[matched]]  
po1$cdd[matched]<- clim1$Cdd[match1[matched]]

##match back to dat
datm= dat[duplicated(dat$siteyear)==FALSE,]
match1= match(po1$siteyear, datm$siteyear)
matched= which(!is.na(match1))
po1$Tmean[matched]<- datm$Tmean[match1[matched]]  
po1$cdd[matched]<- datm$cdd_seas[match1[matched]] 

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

#plot
pdf("PhenOverlap_byYear.pdf", height = 10, width = 10)
ggplot(data=po1[which(po1$metric==1),], aes(x=year, y = value, color=site, shape=period))+geom_point() +facet_grid(sp1~sp2, drop=TRUE)+theme_bw() +geom_smooth(method=lm, se=FALSE) #geom_line
dev.off() 

#--------
#plot by temp and gdd

#add elevation
po1$elevation= as.factor(elevs[match(po1$site, sites)])

#drop Chataqua for now 
#po1= po1[which(po1$site!="CHA"),]

#focus on sites B1 and C1 for now
po1= po1[which(po1$site %in% c("B1", "C1")) ,] #[which(po1$metric==2),] #

#plot

#overlap by temp
pdf("PhenOverlap_byTemp.pdf", height = 10, width = 10)
ggplot(data=po1[which(po1$metric==1),], aes(x=Tmean, y = value, color=site, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()#+geom_smooth(method=lm, se=FALSE)
dev.off()

#*********************************************
#overlap by GDD

#assess significance
#calculate significant regressions
#apply through combinations of species and elevations
po1$elevspec= paste(po1$elevation, po1$sp, sep="")
elevspec= matrix(unique(po1$elevspec))

#---
#EXTRACT COEFFS
p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(po1$value[which(po1$metric==1 & po1$elevation==substr(x,1,4)& po1$sp==substr(x,5,nchar(x)))] ~ po1$cdd[which(po1$metric==1 & po1$elevation==substr(x,1,4)&po1$sp==substr(x,5,nchar(x)) )]) )$coefficients[2,])
#combine
p.mat=as.data.frame(cbind(elevspec,t(p.gdd) ))
#make numeric
p.mat$Estimate= as.numeric(as.character(p.mat$Estimate))
p.mat$`Std. Error`= as.numeric(as.character(p.mat$`Std. Error`))
#CIs
p.mat$upCI= p.mat$Estimate + 1.96*p.mat$`Std. Error`
p.mat$lowCI= p.mat$Estimate - 1.96*p.mat$`Std. Error`
#---

#extract p-values
p.gdd= apply(elevspec,1, FUN=function(x) summary(lm(po1$value[which(po1$metric==1 & po1$elevation==substr(x,1,4)& po1$sp==substr(x,5,nchar(x)))] ~ po1$cdd[which(po1$metric==1 & po1$elevation==substr(x,1,4)&po1$sp==substr(x,5,nchar(x)) )]) )$coefficients[2,4])
#combine
p.mat=as.data.frame(cbind(elevspec,p.gdd))
#add columns for significance
p.mat$sig.gdd<-"nonsignificant"
p.mat$sig.gdd[which(p.gdd<0.05)]="significant"

#add back to matrix
match1= match(po1$elevspec, elevspec)
po1$sig.gdd= factor(p.mat[match1,"sig.gdd"], levels=c("significant","nonsignificant"))
#---

pdf("PhenOverlap_byGDD.pdf", height = 8, width = 10)
ggplot(data=po1[which(po1$metric==1),], aes(x=cdd, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+
  facet_grid(sp1~sp2, drop=TRUE)+theme_bw()+geom_smooth(method="lm", se=FALSE, aes(linetype=sig.gdd))+ylab("Phenological overlap")+xlab("Growing degree days")+ 
  scale_color_manual(values=c("darkorange", "blue","darkgreen","purple"))#+xlim(200,750)
dev.off()
#**********************************************

#overlap by GDD by site
pdf("PhenOverlap_byGDD_site.pdf", height = 8, width = 10)
ggplot(data=po1[which(po1$metric==1),], aes(x=cdd, y = value, color=sp))+geom_point(aes(shape=period, fill=period), size=2)+
  facet_wrap(~site, drop=TRUE, scales="free")+theme_bw()+geom_smooth(method="lm", se=FALSE)+ylab("Phenological overlap")+xlab("Growing degree days")+ 
  scale_color_manual(values=rainbow(10))+ylim(0,1)
dev.off()

#just B1 and C1
pdf("PhenOverlap_byGDD_bysite.pdf", height = 8, width = 10)
ggplot(data=po1[which(po1$metric==1 & po1$site %in% c("A1","B1")),], aes(x=cdd, y = value, color=elevation))+geom_point(aes(shape=period, fill=period), size=2)+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()+geom_smooth(method="lm", se=FALSE)+ylab("Phenological overlap")+xlab("Growing degree days")+ scale_color_manual(values=c("darkorange", "blue","darkgreen","purple"))
dev.off()

#===================================
# STATS GROUP BY SPECIES
mod1= lm(value~ cdd*elevation+sp2+sp1, data=po2)
mod1= lm(value~ cdd*elevation+sp, data=po2)

#by focal species
sp.k=6
po3= subset(po2, po2$sp1==specs[sp.k] | po2$sp2==specs[sp.k])
#switch species order
p.temp= po3$sp1
inds= which(po3$sp2==specs[sp.k])
if(length(inds)>0){
 po3$sp1[inds]=specs[sp.k] 
 po3$sp2[inds]=po3$sp1[inds] 
}

po3$elevation= as.numeric(as.character(po3$elevation))
mod1= lm(value~ cdd*elevation+sp2, data=po3)
summary(mod1)

#--------------------------------
#significant trends by species
#Code from Stuart

finalDat=po2

# Create results storage
modelOutput <- array(NA, dim = c(length(unique(finalDat$sp)), 4, length(sites)))
colnames(modelOutput) <- c("slope", "SE", "tvalue", "pvalue")

coms= unique(finalDat$sp)

for(site in 1:length(sites)){
  for(com in 1:length(unique(finalDat$sp))){
    
    # Subset to 1 site and 1 species combination
    dat <- finalDat[finalDat$site == sites[site] & finalDat$sp == coms[com],]
    
    # Remove year 2009 for B1
    if(sites[site] == "B1"){
      dat <- dat[dat$year != 2009,]
    }
    
    # Only try to create model if data exists
    if(sum(!is.na(dat$value)) > 0){
      # Create linear model
      mod <- lm(value ~ cdd, data = dat)
      
      # Extract coefficients table
      modelOutput[com,,site] <- summary(mod)[["coefficients"]][2,]
    }
  }
}

# Combine all output in a data frame
combIDs <- finalDat[match(sort(unique(finalDat$sp)), finalDat$sp), c("sp", "sp1", "sp2")]
coefs <- rbind(modelOutput[,,1], modelOutput[,,2], modelOutput[,,3], modelOutput[,,4])
IDs <- rbind(combIDs, combIDs, combIDs, combIDs)
modelResults <- cbind(IDs, coefs)
modelResults$site <- rep(sites, each = 15)

# View marginally significant results
View(modelResults[!is.na(modelResults$pvalue) & modelResults$pvalue < 0.1,])

#----

