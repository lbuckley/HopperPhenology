library(reshape2)

#LOAD data
setwd("/Volumes/GoogleDrive/My Drive/AlexanderResurvey/DataForAnalysis/grasshoppers/ForLTER/")

#-----------------
#historic

hj= read.csv("HopperData_NymphalHist.csv")

#-----------------
#juveniles
b1.j= read.csv("HopperData_Sep2019_NymphalB1_juvenile.csv")
c1.j= read.csv("HopperData_Sep2019_NymphalC1_juvenile.csv")
chaut.j= read.csv("HopperData_Sep2019_NymphalChaut_juvenile.csv")
#add sites
b1.j$site="B1"
c1.j$site="C1"
chaut.j$site="Chaut"
#bind
nj= rbind(b1.j, c1.j, chaut.j)
#add date
date= as.Date(nj$date, format="%m/%d/%y")
nj$ordinal= format(date,"%j")
nj$year= format(date,"%Y")
nj$sp.site.yr.ord= paste(nj$Species, nj$site, nj$year, nj$ordinal, sep="_")

#-----------------
#adults
b1.a= t(read.csv("HopperData_Sep2019_NymphalB1_adult.csv", header=FALSE) )
c1.a= t(read.csv("HopperData_Sep2019_NymphalC1_adult.csv", header=FALSE) )
chaut.a= t(read.csv("HopperData_Sep2019_NymphalChaut_adult.csv", header=FALSE) )
#fix format
colnames(b1.a)<- b1.a[1,]
colnames(c1.a)<- c1.a[1,]
colnames(chaut.a)<- chaut.a[1,]
b1.a=b1.a[-1,]
c1.a=c1.a[-1,]
chaut.a=chaut.a[-1,]
b1.a= as.data.frame(b1.a)
c1.a= as.data.frame(c1.a)
chaut.a= as.data.frame(chaut.a)
#to long format
b1.a.m= melt(b1.a, id = c("DATE","Ordinal Date:"))
c1.a.m= melt(c1.a, id = c("DATE","OrdinalDate"))
chaut.a.m= melt(chaut.a, id = c("Date","Ordinal Date:"))
colnames(b1.a.m)<- colnames(c1.a.m)
colnames(chaut.a.m)<- colnames(c1.a.m)
#add sites
b1.a.m$site="B1"
c1.a.m$site="C1"
chaut.a.m$site="Chaut"
#bind
nad= rbind(b1.a.m, c1.a.m, chaut.a.m)

#add date
date= as.Date(nad$DATE, format="%m/%d/%y")
nad$ordinal= format(date,"%j")
nad$year= format(date,"%Y")
nad$sp.site.yr.ord= paste(nad$variable, nad$site, nad$year, nad$ordinal, sep="_")

#match format
colnames(nj)[6]<-"X6.Adult"
nad$GDDs<-NA
nad1= cbind(nad[,c("year","DATE","ordinal","GDDs","variable","value")], matrix(NA, nrow(nad),5), nad[,c("site","sp.site.yr.ord")] )   
colnames(nad1)<- colnames(nj)

#--------------
#match
nad1$X6.Adult= as.numeric(as.character(nad1$X6.Adult))
nj$X6.Adult= NA
nj$X6.Adult= as.numeric(nj$X6.Adult)

match1=match(nj$sp.site.yr.ord, nad1$sp.site.yr.ord)
matched= which(!is.na(match1))

nj$X6.Adult[matched]= nad1$X6.Adult[match1[matched]] 

#add unmatched
match2=match(nad1$sp.site.yr.ord, nj$sp.site.yr.ord)
unmatched= which(is.na(match2))

nall= rbind(nj, nad1[unmatched,])

#--------------
#WRITE OUT

write.csv(nall, "NymphalData.csv")

#++++++++++++++++++++=

fdir= "/Volumes/GoogleDrive/My\ Drive/AlexanderResurvey/DataForAnalysis/"

#load climate data
setwd( paste(fdir, "climate", sep="") )   
clim= read.csv("AlexanderClimateAll_filled_May2018.csv")

#load hopper data
setwd( paste(fdir, "grasshoppers/SexCombined/", sep="") )
hop= read.csv("HopperData_May2018.csv")

#----------
#Add climate data to nall
nall$sjy= paste(nall$site, nall$ordinal, nall$year, sep="_" )

match1= match(nall$sjy, hop$sjy)
