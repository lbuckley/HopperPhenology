library(reshape)
library(tidyr)

#data currently from phenology file
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

#----------------------------
for(year in 1:length(years)){
  for(sp in 1:length(specs)){
    
    #sum all individuals by species and year 
    dat.sub= dat[dat$year==years[year],]
    dat.sub= dat.sub[ which(dat.sub$species==specs[sp]),]
    #restrict to adults
    dat.sub= subset(dat.sub, dat.sub$in6>0 )
    
    if( nrow(dat.sub)>0 ){ #check data exists
      
      ords= unique(dat.sub$ordinal)
      
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

#---------------------------------------------
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

#----------------------------------------------
#Plot days overlab by elevation and species

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
