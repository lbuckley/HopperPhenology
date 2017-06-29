#load libraries
library(ggplot2)
library(plyr) 
library(dplyr)
library(reshape)
library(tidyr)

sites= c("CHA", "A1", "B1", "C1", "D1")  #Redfox: 1574
elevs= c(1752, 2195, 2591, 3048, 3739)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

#======================================================
#READ DATA

fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#---------------------
#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd_sum = cumsum(dd_sum),cdd_june = cumsum(dd_june),cdd_july = cumsum(dd_july),cdd_aug = cumsum(dd_aug),cdd_early = cumsum(dd_early),cdd_mid = cumsum(dd_mid),cdd_ac = cumsum(dd_ac),cdd_mb = cumsum(dd_mb),cdd_ms = cumsum(dd_ms) ) 

#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )
hop= read.csv("HopperData.csv")

#------------------------------------
#Phenological overlap across elevation

#subset to focal species
specs= c("Aeropedellus clavatus","Camnula pellucida","Melanoplus boulderensis","Melanoplus sanguinipes")

years= unique(na.omit(hop$year))
sites= unique(hop$site)

dat=hop

#all sites combinations
sp.combs= combn(sites, 2)

#set up data storage
po= array(data = NA, dim = c(length(sp.combs[1,]), length(years), length(specs)) )

#combine species combinations names
sp.combs.n= paste(sp.combs[1,],sp.combs[2,],sep="_")

#add po names
dimnames(po)[[1]] <- sp.combs.n
dimnames(po)[[2]] <- years
dimnames(po)[[3]] <- specs

#group ordinal dates
getMondays <- function(year) {
  days <- as.POSIXlt(paste(year, 70:300, sep="-"), format="%Y-%j")
  Ms <- days$yday[days$wday==1]
  Ms[!is.na(Ms)]  # Needed to remove NA from day 366 in non-leap years
}


#----------------------------
for(year in 1:length(years)){
  
  dat.sub.all= dat[dat$year==years[year],]
  #group ordinals
  doys= getMondays(years[year])
  dat.sub.all$week= cut(dat.sub.all$ordinal, doys, labels=FALSE )
  
  for(sp in 1:length(specs)){
    
    #sum all individuals by species and year 
    dat.sub= dat.sub.all[ which(dat.sub.all$species==specs[sp]),]
    #restrict to adults
    dat.sub= subset(dat.sub, dat.sub$in6>0 )
    
    if( nrow(dat.sub)>0 ){ #check data exists
      
      ords= unique(dat.sub$week) #dat.sub$ordinal
      
      #totals by ordinal date
      dat.agg= dat.sub[,c("site", "in6")] #extract adults
      dat.agg$site= as.factor(dat.sub$site)
      dat.agg= aggregate(dat.agg$in6, by=list(dat.agg$site), FUN=sum) #sum number adults
      names(dat.agg)= c("site", "YrTot")
      
      #turn into proportions
      match1= match( dat.sub$site, dat.agg$site)
      dat.sub$YrTot= dat.agg$YrTot[match1] #Add total number adults per species per year
      dat.sub= dat.sub[which(dat.sub$YrTot>0),]
      dat.sub$AdultProp= dat.sub$in6/dat.sub$YrTot #Add proportion of annual individuals for each ordinal date, normalize species with different abundances 
      #subset to dates with adults
      dat.sub= dat.sub[which(dat.sub$AdultProp>0),]
      
      #find species combinations potentially present
      sites1= unique(dat.sub$site)
      #find species combinations present
      inds=which(sp.combs[1,]%in%sites1 & sp.combs[2,]%in%sites1)
      
      #loop through site combinations
      if(length(inds)>1) for(ind.k in 1:length(inds) ){
        sp1= sp.combs[1,inds[ind.k] ]
        sp2= sp.combs[2,inds[ind.k] ]
        #subset to combinations
        dat.sub2= dat.sub[which(dat.sub$site==sp1 | dat.sub$site==sp2),]
        ords1= unique(dat.sub2$ordinal) #ordinal dates with adults present
        
        num=0
        dem_pij=0
        dem_pik=0
        
        #loop through ordinal dates
        for(j in 1:length(ords1) ){
          dat.ord= subset(dat.sub2, dat.sub2$ordinal==ords1[j] ) #subset to ordinal date
          prop1= dat.ord[match(sp1, dat.ord$site),"AdultProp"]
          prop2= dat.ord[match(sp2, dat.ord$site),"AdultProp"]
          #replace NAs with zero
          if(is.na(prop1))prop1=0
          if(is.na(prop2))prop2=0
          
          #new metric?: sum across dates: (max proportion-overlap)/total, On the Analysis of Phenological Overlap
          #po[inds[c], year, site]= po[inds[c], year, site]+ (max(prop1,prop2)-abs(prop1-prop2))/(prop1+prop2)
        
          #Pianka (1974) metric, http://www.jstor.org/stable/pdf/4217327.pdf, 
          num= num+prop1*prop2
          dem_pij= dem_pij + prop1^2
          dem_pik= dem_pik + prop2^2
          
          } #end loop dates
        
        #calc metric
        po[inds[ind.k], year, sp]= num/sqrt(dem_pij*dem_pik)
        
        
      } #end loop species combinations
      
      
    } #end check data
  } #end loop year
} #end loop sites

#-----------------------------------
#PLOTS
# PO: sp comb x years x sites

po1= melt(po, varnames=c("site","year","sp"))

#add time period
po1$period="initial"
po1[which(po1$year>1960),"period"]<-"resurvey"

#Change years for plotting ### FIX
po1[which(po1$year==1958),"year"]= 1958+40
po1[which(po1$year==1959),"year"]= 1959+40
po1[which(po1$year==1960),"year"]= 1960+40

#separate species names for plotting
po2 = transform(po1, site = colsplit(site, split = "_", names = c('1', '2')))
po1$sp1=po2$site$X1
po1$sp2=po2$site$X2

#make site factor for plotting
po1$sp1= factor(po1$sp1, levels= c("CHA","A1","B1","C1") )
po1$sp2= factor(po1$sp2, levels= c("CHA","A1","B1","C1") )

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")
pdf("PopPhenOverlap_byYear.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=year, y = value, color=sp ))+geom_point() +facet_grid(sp1~sp2, drop=TRUE)+theme_bw() #+geom_smooth(method=lm) #geom_line
dev.off() 

#--------------------------------
#PLOT LINES CONNECTING YEARS

#add elevation
po1$elev1=elevs[match(po1$sp1,sites)]
po1$elev2=elevs[match(po1$sp2,sites)]
po1$site_year_sp= paste(po1$site, po1$year,po1$sp, sep="")

#restrict to recent years
po2= po1[which(po1$year>2005),]

#combine first and last
po3=po2[,c(1:4,8:10)]

#drop NAs
po3= po3[which(!is.na(po3$value)),]

po3= po3 %>% gather(key = elevs, value = elevation, -site, -year,-sp, -value,-site_year_sp)
po3$year= as.factor(po3$year)

#### USE
phenoverlap.plot=ggplot(data=po3, aes(x=elevation, y = value, color=sp, group=site_year_sp))+geom_line(size=1)+theme_bw()+ylab("overlap") #linetype=year, 

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")
pdf("PopPhenOverlap_Lines_NoLinetype.pdf", height = 6, width = 6)
ggplot(data=po3, aes(x=elevation, y = value, color=sp, group=site_year_sp))+geom_line(size=1)+theme_bw()+ylab("overlap") #linetype=year, 
dev.off()
#--------------------------------

#Average across years
#restrict to recent years
po2= po1[which(po1$year>2005),]
po2= po2 %>% group_by(sp,site,sp, sp1, sp2) %>% summarise(N= length(na.omit(value)), value=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE), se = sd / sqrt(N))
po2$value[is.nan(po2$value)] = NA

ggplot(data=po2, aes(x=sp, y = value, color=sp, group=site_year_sp))+geom_point()+geom_line()+theme_bw() #+geom_smooth(method=lm) #geom_line

po3= po2[which(po2$sp1=="A1"),]
#add elevation
po2$elev=elevs[match(po2$sp2,sites)]

plot1=ggplot(data=po2, aes(x=elev, y = value, color=sp))+geom_point(aes(shape=sp1), size=3)+geom_smooth(method=lm, se=FALSE) +theme_bw()+ theme(legend.position="right")

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")
pdf("PopPhenOverlap_YearAve.pdf", height = 6, width = 6)
plot1
dev.off() 

#--------
#plot by temp and gdd
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
po1$siteyear1= paste(po1$sp1, po1$year, sep="")
po1$siteyear2= paste(po1$sp2, po1$year, sep="")

po1$Tmean1=NA
po1$cdd1=NA
po1$Tmean2=NA
po1$cdd2=NA

#select cdd metric
clim1=clim
clim1$Cdd=  clim1$Cdd_july

match1= match(po1$siteyear1, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean1[matched]<- clim1$Mean[match1[matched]]  
po1$cdd1[matched]<- clim1$Cdd[match1[matched]]
match1= match(po1$siteyear2, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean2[matched]<- clim1$Mean[match1[matched]]  
po1$cdd2[matched]<- clim1$Cdd[match1[matched]]
#Average site values
po1$Tmean=rowMeans(po1[,c("Tmean1", "Tmean2")])
po1$cdd=rowMeans(po1[,c("cdd1", "cdd2")])

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#overlap by temp
pdf("PopPhenOverlap_byTemp.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=Tmean, y = value, color=sp, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()#+geom_smooth(method=lm, se=FALSE)
dev.off()

#overlap by GDD
pdf("PopPhenOverlap_byGDD.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=cdd, y = value, color=sp, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()
dev.off()

#===========================================
#Analyze by species

#check first and last sampling dates across sites
seas.samp= ddply(hop1, .(site,year), summarize, first=min(ordinal,na.rm=TRUE), last=max(ordinal,na.rm=TRUE) )

#Plot Season by species
#Elev for A1  B1  C1  CHA
elevs= c(2195, 2591, 3048, 1752)

#restrict to dates with adults
hop1= hop[which(hop$in6>0),]

#Find first and last adult date by species, site, year
seas= ddply(hop1, .(species,site,year), summarize, first=min(ordinal,na.rm=TRUE), last=max(ordinal,na.rm=TRUE) )

#restrict to focal species
seas= seas[which(seas$species %in% specs),]
seas$ind= 1:nrow(seas)

#combine first and last
seas1= seas %>% gather(key = per, value = ordinal, -species, -site, -year, -ind)

#add elevation
match1= match(seas1$site, sites)
seas1$elev= elevs[match1]

#jitter years
years1= sort(years)
seas1$elev.adj= seas1$elev + match(seas1$year, years1)*10

#season by species
pdf("SeasByElev.pdf", height = 10, width = 10)
ggplot(data=seas1, aes(x=ordinal, y = elev.adj, color=year, group=ind))+geom_line()+facet_grid(.~species)+theme_bw()#+geom_smooth(method=lm, se=FALSE) 
dev.off()
#jitter by year

#-------------
## calculate median across individuals
hop1= hop1[order(hop1$ordinal),]

#restrict to focal species
hop1= hop1[which(hop1$species %in% specs),]

hop1= hop1[which(hop1$year>2005),]

#cumulative sum of individuals within groups
hop1 = hop1 %>% group_by(species,site,year) %>% arrange(species,site,year,ordinal) %>% mutate(csind = cumsum(in6))

#number of median individual
hop3 = hop1 %>% group_by(species,site,year) %>% arrange(species,site,year,ordinal) %>% mutate(medind = max(csind)/2, q20ind=max(csind)*0.2, q80ind=max(csind)*0.8)

#date of median individual
hop3$inddif= abs(hop3$medind-hop3$csind) #difference from median individual
hop3$inddif.q20= abs(hop3$q20ind-hop3$csind) #difference from q20 individual
hop3$inddif.q80= abs(hop3$q80ind-hop3$csind) #difference from q80 individual

hop4= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif),]))
hop4.q20= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif.q20),]))
hop4.q80= do.call(rbind,lapply(split(hop3,list(hop3$species, hop3$site, hop3$year)),function(chunk) chunk[which.min(chunk$inddif.q80),]))

#plot
ggplot(data=hop4, aes(x=year, y = ordinal, color=site ))+geom_point()+geom_line()+facet_wrap(~species, ncol=2) +theme_bw()

#-----------
#add elevation
match1= match(hop4$site, sites)
hop4$elev= elevs[match1]
match1= match(hop4.q20$site, sites)
hop4.q20$elev= elevs[match1]
match1= match(hop4.q80$site, sites)
hop4.q80$elev= elevs[match1]

medord_byelev=ggplot(data=hop4, aes(x=elev, y = ordinal, color=year, group=year))+facet_wrap(~species, ncol=2)+theme_bw()+ylab("median doy")+geom_line(lwd=1)#+geom_smooth(method="lm", se=FALSE)

#plot 20th and 80th percentiles
ggplot(data=hop.qs1, aes(x=elev, y = ordinal, color=year, shape=q,group=year_sp))+facet_wrap(~species, ncol=2)+theme_bw()+ylab("20th, 80th doy")+geom_line(lwd=1)+geom_point()

ggplot(data=hop4.q80, aes(x=elev, y = ordinal, color=year, group=year))+facet_wrap(~species, ncol=2)+theme_bw()+ylab("median doy")+geom_line(lwd=1)
ggplot(data=hop4.q20, aes(x=elev, y = ordinal, color=year, group=year))+facet_wrap(~species, ncol=2)+theme_bw()+ylab("median doy")+geom_line(lwd=1)

#combine
hop4.q80$q=80
hop4.q20$q=20
hop5= rbind(hop4.q80, hop4.q20)
hop5$q= as.factor(hop5$q)
hop5$yr_q= paste(hop.qs1$year, hop.qs1$q, sep="_")
hop5$year= as.factor(hop5$year)


#### USE
doyq.plot= ggplot(data=hop5, aes(x=elev, y = ordinal, color=year, linetype=q, group=yr_q))+facet_wrap(~species, nrow=2)+theme_bw()+ylab("doy")+geom_line(lwd=1)

#PLOT median ordinal date by elevation
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperGenetics\\figures\\")
pdf("DOY_byElev.pdf", height = 10, width = 10)
doyq.plot
dev.off()


#---------------------
#First & last
seas= ddply(hop1, .(species,site,year), summarize, first=min(ordinal,na.rm=TRUE), last=max(ordinal,na.rm=TRUE) ) #, q10=quantile(ordinal, probs = 0.1),  q90=quantile(ordinal, probs = 0.9)
#add elevation
match1= match(seas$site, sites)
seas$elev= elevs[match1]

ggplot(data=seas, aes(x=elev, y = first, color=year, group=year))+facet_wrap(~species, ncol=2)+theme_bw()+ylab("first doy")+geom_line(lwd=1)
ggplot(data=seas, aes(x=elev, y = last, color=year, group=year))+facet_wrap(~species, ncol=2)+theme_bw()+ylab("first doy")+geom_line(lwd=1)

#calculate difference in median
hop5= hop4[,c("ordinal","species","year","site")]

#combine first and last
#hop5$site_spec= paste(hop5$site, hop5$species, sep="_")
hop6= hop5 %>% spread(key = site, value = ordinal)
hop6$C1_B1


#----------------------------------------------
#Plot days overlap by elevation and species

#calculate duration of overlap by site combination

seas$sp.si= paste(seas$species,seas$site,sep="_")

seas %>% group_by(species, year) %>% summarise( min   )

#=============================

#subset to focal species
specs= c("Aeropedellus clavatus","Camnula pellucida","Melanoplus boulderensis","Melanoplus sanguinipes")

years= unique(na.omit(seas$year))
sites= unique(seas$site)

dat=seas

#all sites combinations
sp.combs= combn(sites, 2)

#set up data storage
po= array(data = NA, dim = c(length(sp.combs[1,]), length(years), length(specs)) )

#combine species combinations names
sp.combs.n= paste(sp.combs[1,],sp.combs[2,],sep="_")

#add po names
dimnames(po)[[1]] <- sp.combs.n
dimnames(po)[[2]] <- years
dimnames(po)[[3]] <- specs

#----------------------------
for(year in 1:length(years)){
  for(sp in 1:length(specs)){
    
    #sum all individuals by species and year 
    dat.sub= dat[dat$year==years[year],]
    dat.sub= dat.sub[ which(dat.sub$species==specs[sp]),]
   
    if( nrow(dat.sub)>0 ){ #check data exists
      
      #find species combinations potentially present
      sites1= unique(dat.sub$site)
      #find species combinations present
      inds=which(sp.combs[1,]%in%sites1 & sp.combs[2,]%in%sites1)
      
      #loop through site combinations
      if(length(inds)>1) for(ind.k in 1:length(inds) ){
        sp1= sp.combs[1,inds[ind.k] ]
        sp2= sp.combs[2,inds[ind.k] ]
        #subset to combinations
        
        #calculate duration
        po[inds[ind.k], year, sp]= min(dat.sub[which(dat.sub$site==sp1),"last"], dat.sub[which(dat.sub$site==sp2),"last"])-max(dat.sub[which(dat.sub$site==sp1),"first"], dat.sub[which(dat.sub$site==sp2),"first"])
        
      } #end loop species combinations
      
      
    } #end check data
  } #end loop year
} #end loop sites

#-----------------------------------
#PLOTS
# PO: sp comb x years x sites

po1= melt(po, varnames=c("site","year","sp"))

#add ti me period
po1$period="initial"
po1[which(po1$year>1960),"period"]<-"resurvey"

#Change years for plotting ### FIX
po1[which(po1$year==1958),"year"]= 1958+40
po1[which(po1$year==1959),"year"]= 1959+40
po1[which(po1$year==1960),"year"]= 1960+40

#separate species names for plotting
po2 = transform(po1, site = colsplit(site, split = "_", names = c('1', '2')))
po1$sp1=po2$site$X1
po1$sp2=po2$site$X2

#make site factor for plotting
po1$sp1= factor(po1$sp1, levels= c("CHA","A1","B1","C1") )
po1$sp2= factor(po1$sp2, levels= c("CHA","A1","B1","C1") )

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")
pdf("DurOverlap.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=year, y = value, color=sp ))+geom_point() +facet_grid(sp1~sp2, drop=TRUE)+theme_bw() #+geom_smooth(method=lm) #geom_line
dev.off() 

#--------
#plot by temp and gdd
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
po1$siteyear1= paste(po1$sp1, po1$year, sep="")
po1$siteyear2= paste(po1$sp2, po1$year, sep="")

po1$Tmean1=NA
po1$cdd1=NA
po1$Tmean2=NA
po1$cdd2=NA
match1= match(po1$siteyear1, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean1[matched]<- clim1$Mean[match1[matched]]  
po1$cdd1[matched]<- clim1$Cdd[match1[matched]]
match1= match(po1$siteyear2, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean2[matched]<- clim1$Mean[match1[matched]]  
po1$cdd2[matched]<- clim1$Cdd[match1[matched]]
#Average site values
po1$Tmean=rowMeans(po1[,c("Tmean1", "Tmean2")])
po1$cdd=rowMeans(po1[,c("cdd1", "cdd2")])

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#overlap by temp
pdf("DurOverlap_byTemp.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=Tmean, y = value, color=sp, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()#+geom_smooth(method=lm, se=FALSE)
dev.off()

#overlap by GDD
pdf("DurPhenOverlap_byGDD.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=cdd, y = value, color=sp, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()
dev.off()

#---------------------------------------------
