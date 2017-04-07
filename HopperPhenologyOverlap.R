library(reshape)
library(tidyr)

#data currently from phenology file
#------------------------------------

#subset to focal species
specs= c("Aeropedellus clavatus","Chloealtis abdominalis","Camnula pellucida","Melanoplus dawsoni","Melanoplus boulderensis","Melanoplus sanguinipes")

years= unique(na.omit(hop$year))
sites= unique(hop$site)

dat=hop

#all species combinations
sp.combs= combn(specs, 2)

#set up data storage
po= array(data = 0, dim = c(length(sp.combs[1,]), length(years), length(sites)) )

#combine species combinations names
sp.combs.n= paste(sp.combs[1,],sp.combs[2,],sep="_")

#add po names
dimnames(po)[[1]] <- sp.combs.n
dimnames(po)[[2]] <- years
dimnames(po)[[3]] <- sites

#----------------------------
for(year in 1:length(years)){
  for(site in 1:length(sites)){
    
    #sum all individuals by species and year 
    dat.sub= dat[dat$year==years[year],]
    dat.sub= dat.sub[ which(dat.sub$site==sites[site]),]
    #restrict to adults
    dat.sub= subset(dat.sub, dat.sub$in6>0 )
    
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
        ords1= unique(dat.sub2$ordinal) #ordinal dates with adults present
        
        num=0
        dem_pij=0
        dem_pik=0
        
        #loop through ordinal dates
        for(j in 1:length(ords1) ){
          dat.ord= subset(dat.sub2, dat.sub2$ordinal==ords1[j] ) #subset to ordinal date
          prop1= dat.ord[match(sp1, dat.ord$species),"AdultProp"]
          prop2= dat.ord[match(sp2, dat.ord$species),"AdultProp"]
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
        po[inds[ind.k], year, site]= num/sqrt(dem_pij*dem_pik)
        
        
      } #end loop species combinations
      
      
    } #end check data
  } #end loop year
} #end loop sites

#-----------------------------------
#PLOTS
# PO: sp comb x years x sites

po1= melt(po, varnames=c("sp","year","site"))

#add time period
po1$period="initial"
po1[which(po1$year>1960),"period"]<-"resurvey"

#Change years for plotting ### FIX
po1[which(po1$year==1958),"year"]= 1958+40
po1[which(po1$year==1959),"year"]= 1959+40
po1[which(po1$year==1960),"year"]= 1960+40

#separate species names for plotting
po2 = transform(po1, sp = colsplit(sp, split = "_", names = c('1', '2')))
po1$sp1=po2$sp$X1
po1$sp2=po2$sp$X2

#order by early and late species
#ave phen
#hop.el = hop1 %>% group_by(species,site) %>% arrange(species,site) %>% mutate(phen = mean(ordinal))
hop.agg= aggregate(hop1, list(hop1$species),FUN=mean) #fix for hop1$site,
hop.agg= hop.agg[order(hop.agg$ordinal),]

#make species factor for plotting
po1$sp1= factor(po1$sp1, levels=hop.agg$Group.1)
po1$sp2= factor(po1$sp2, levels=hop.agg$Group.1)

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")
pdf("PhenOverlap_byYear.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=year, y = value, color=site ))+geom_point() +facet_grid(sp1~sp2, drop=TRUE)+theme_bw() #+geom_smooth(method=lm) #geom_line
dev.off() 

#--------
#plot by temp and gdd
clim1$siteyear= paste(clim1$Site, clim1$Year, sep="")
po1$siteyear= paste(po1$site, po1$year, sep="")

po1$Tmean=NA
po1$cdd=NA
match1= match(po1$siteyear, clim1$siteyear)
matched= which(!is.na(match1))
po1$Tmean[matched]<- clim1$Mean[match1[matched]]  
po1$cdd[matched]<- clim1$Cdd[match1[matched]]

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#overlap by temp
pdf("PhenOverlap_byTemp.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=Tmean, y = value, color=site, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()#+geom_smooth(method=lm, se=FALSE)
dev.off()

#overlap by GDD
pdf("PhenOverlap_byGDD.pdf", height = 10, width = 10)
ggplot(data=po1, aes(x=cdd, y = value, color=site, shape=period))+geom_point()+facet_grid(sp1~sp2, drop=TRUE)+theme_bw()
dev.off()

