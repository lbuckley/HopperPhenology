#load libraries
library(ggplot2)
library(plyr) 
library(dplyr)

sites= c("CHA", "A1", "B1", "C1", "D1")  #Redfox: 1574
elevs= c(1752, 2195, 2591, 3048, 3739)

#source degree days function
setwd("C:\\Users\\Buckley\\Documents\\HopperPhenology\\")
source("degreedays.R")

#--------------------------------------
fdir= "C:\\Users\\Buckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"
#fdir= "C:\\Users\\lbuckley\\Google Drive\\AlexanderResurvey\\DataForAnalysis\\"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled.csv")

#---------------------
#cummulative degree days
#cumsum within groups
clim = clim %>% group_by(Year,Site) %>% arrange(Julian) %>% mutate(cdd_sum = cumsum(dd_sum),cdd_june = cumsum(dd_june),cdd_july = cumsum(dd_july),cdd_aug = cumsum(dd_aug),cdd_early = cumsum(dd_early),cdd_mid = cumsum(dd_mid),cdd_ac = cumsum(dd_ac),cdd_mb = cumsum(dd_mb),cdd_ms = cumsum(dd_ms) ) 

#========================================================================
#load hopper data
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )

hop.b1= read.csv( "B1_1959-2015_eggDiapause.csv" )
hop.c1= read.csv( "C1_1959-2015_eggDiapause.csv" )
hop.cha= read.csv( "CHA_1959-2012_eggDiapause.csv" )
#add old a1 data
hop.a1= read.csv( "A1data_old.csv" )

#update A1 data
hop.a1.up= read.csv( "A1_1958-2011_allSpecies.csv" ) 
#combine for now across sex and subsite
hop.a1.up2= ddply(hop.a1.up, .(species,year,ordinal), summarize, date=date[1], GDDs= GDDs[1], in6= sum(in6),in5= sum(in5),in4= sum(in4),in3= sum(in3),in2= sum(in2),in1= sum(in1),total= sum(total), month=month[1], day=day[1], site=site[1] )
hop.a1= hop.a1.up2[,c("date","ordinal","GDDs","species","in6","in5","in4","in3","in2","in1","total","month","day","year","site")]                  

#add site
hop.b1$site="B1"
hop.c1$site="C1"
hop.cha$site="CHA"

#combine
hop= rbind(hop.a1,hop.b1,hop.c1,hop.cha)

#add time period
hop$period="initial"
hop[which(hop$year>1960) ,"period"]="resurvey"

#--------------------------------------
#CHECK DATA

#check years
sort(unique(clim$Year)) #through 2014
sort(unique(hop$year)) #through 2015

#remove white space in species names
hop$species=trimws(hop$species, which = "both")
#check species
sort(unique(hop$species))
#Fix typos
hop[which(hop$species=="Melanoplus bouderensis"),"species"]="Melanoplus boulderensis"
hop[which(hop$species=="Melanoplus dodgei"),"species"]="Melanoplus boulderensis"
hop[which(hop$species=="Melanoplus bivitattus"),"species"]="Melanoplus bivittatus"
hop[which(hop$species=="Cratypedes neglectus"),"species"]="Cratypledes neglectus"
#replace _ with space
hop$species=gsub("_"," ", hop$species)

#------------------------------------
#Add GDD data

#columns for matching
clim$sjy= paste(clim$Site, clim$Julian, clim$Year, sep="_")
hop$sjy= paste(hop$site, hop$ordinal, hop$year, sep="_")

#match
match1= match(hop$sjy,clim$sjy)
matched= which(!is.na(match1))

hop$dd_sum=NA;hop$cdd_sum=NA; hop$cdd_june=NA; hop$cdd_july=NA; hop$cdd_aug=NA; hop$cdd_early=NA; hop$cdd_mid=NA; hop$cdd_ac=NA; hop$cdd_mb=NA; hop$cdd_ms=NA

hop[matched,c("dd_sum","cdd_sum","cdd_june","cdd_july","cdd_aug","cdd_early","cdd_mid","cdd_ac","cdd_mb","cdd_ms")]= clim[match1[matched], c("dd_sum","cdd_sum","cdd_june","cdd_july","cdd_aug","cdd_early","cdd_mid","cdd_ac","cdd_mb","cdd_ms")]
#-------------------------------------

#Write hopper data out
setwd( paste(fdir, "grasshoppers\\SexCombined\\", sep="") )
#write.csv(hop, "HopperData.csv")

#======================================================
#ASSESS A1 data

#combine for now across sex and subsite
hop.a1.site= ddply(hop.a1.up, .(species,year,ordinal, subsite), summarize, date=date[1], GDDs= GDDs[1], in6= sum(in6),in5= sum(in5),in4= sum(in4),in3= sum(in3),in2= sum(in2),in1= sum(in1),total= sum(total), month=month[1], day=day[1], site=site[1] )
#hop.a1.site= hop.a1.site[,c("date","ordinal","GDDs","species","in6","in5","in4","in3","in2","in1","total","month","day","year","site")] 

#-------------------------------------
#DI

dat=hop.a1.site

#drop unknown
dat= dat[-which(dat$subsite=="Unknown"),]

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

#replace _
dat$species= gsub("_"," ", dat$species)

#reduce to focal species
dat= dat[dat$species %in% specs,]

#code period
dat$per=1
dat$per[dat$year>1960]=2

#------------------------------------------------------------
#PLOTS

dat$species= as.factor(dat$species)
dat$year= as.factor(dat$year)
dat$per= as.factor(dat$per)

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#--------------------
#DEVELOPMENTAL INDEX
#Plot DI by ordinal date
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=year, linetype=subsite)) +geom_point(aes(shape=per), size=2)+theme_bw()+geom_line() +facet_grid(species~.)

#Plot DI by GDD
dat$cdd= dat$cdd_sum
di.plot= ggplot(data=dat, aes(x=cdd, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+geom_line()+xlim(0,600)
#note xlim restricted

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#min ordinal date by dd, #plot regression
pdf("DI_A1.pdf", height = 10, width = 5)
di.plot
dev.off()

#no sig difference by subsite
mod1= lm(dat$DI~ dat$species + dat$year * dat$ordinal * dat$subsite )

#------------------
#PHENOLOGY PLOT

#subset to dates with adults
dat1= dat[which(dat$in6>0),]

dat1$year= as.numeric(as.character(dat1$year))
dat1= dat1[which(dat1$species %in% specs ),]
dat1$period= "initial"
dat1$period[which(dat1$year>1960)]= "resurvey"

#metrics across years
hop2= ddply(dat1, c("site","subsite", "year","species","period"), summarise,
            min = min(ordinal, na.rm=TRUE), mean = mean(ordinal, na.rm=TRUE) )

#min ordinal date
ggplot(data=hop2, aes(x=year, y = min, color=subsite, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date
ggplot(data=hop2, aes(x=year, y = mean, color=subsite, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()

#====================================================
# CHECK SEXES

#combine for now across sex and subsite
hop.a1.site= ddply(hop.a1.up, .(species,year,ordinal, gender), summarize, date=date[1], GDDs= GDDs[1], in6= sum(in6),in5= sum(in5),in4= sum(in4),in3= sum(in3),in2= sum(in2),in1= sum(in1),total= sum(total), month=month[1], day=day[1], site=site[1] )

#-------------------------------------
#DI

dat=hop.a1.site

#drop unknown sex
dat= dat[-which(dat$gender %in% "unknown" ),]

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

#replace _
dat$species= gsub("_"," ", dat$species)

#reduce to focal species
dat= dat[dat$species %in% specs,]

#code period
dat$per=1
dat$per[dat$year>1960]=2

#------------------------------------------------------------
#PLOTS

dat$species= as.factor(dat$species)
dat$year= as.factor(dat$year)
dat$per= as.factor(dat$per)

#get rid of entries with low sample size: at least three individuals, #but keep if all first instar or adult
drop= which(dat$total<3) # & dat$DI!=1 & dat$DI!=6
if(length(drop)>0) dat=dat[-drop,]

#--------------------
#DEVELOPMENTAL INDEX
#Plot DI by ordinal date
di.plot= ggplot(data=dat, aes(x=ordinal, y = DI, color=year, linetype=gender)) +geom_point(aes(shape=per), size=2)+theme_bw()+geom_line() +facet_grid(species~.)

#Plot DI by GDD
dat$cdd= dat$cdd_sum
di.plot= ggplot(data=dat, aes(x=cdd, y = DI, color=year))+facet_grid(species~elev) +geom_point(aes(shape=period), size=2)+theme_bw()+geom_line()+xlim(0,600)
#note xlim restricted

#plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\GrasshopperPhenSynch\\figures\\")

#min ordinal date by dd, #plot regression
pdf("DI_A1_bySex.pdf", height = 10, width = 5)
di.plot
dev.off()

#no sig difference by sex
mod1= lm(dat$DI~ dat$species + dat$year * dat$ordinal * dat$gender )

#------------------
#PHENOLOGY PLOT

#subset to dates with adults
dat1= dat[which(dat$in6>0),]

dat1$year= as.numeric(as.character(dat1$year))
dat1= dat1[which(dat1$species %in% specs ),]
dat1$period= "initial"
dat1$period[which(dat1$year>1960)]= "resurvey"

#metrics across years
hop2= ddply(dat1, c("site","gender", "year","species","period"), summarise,
            min = min(ordinal, na.rm=TRUE), mean = mean(ordinal, na.rm=TRUE) )

#min ordinal date
ggplot(data=hop2, aes(x=year, y = min, color=gender, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()
#mean ordinal date
ggplot(data=hop2, aes(x=year, y = mean, color=gender, shape=period))+geom_point()+geom_line()+facet_wrap(~species, ncol=3) +theme_bw()




